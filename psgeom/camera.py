#!/usr/bin/env python

"""
Detector Geometries
===================

Briefly, geometries are organized heriarchically: there are `SensorElements`
that compose the leaves of a tree. These represent arbitrary objects that
actually measure e.g. photon intensities of scattered x-rays. An example of
such an element is an ASIC on a Pilatus, or a 2x1 on a CSPAD.

These elements are then "mounted" onto `CompoundCamera`s, which represent
physical units that translate and rotate together. For example, 8 CSPAD 2x1s
are mounted on a single "Quad", that moves as a collective unit. The elements
composing a CompoundCamera` instance can be SensorElements or other 
CompoundCameras.

A note about how the heirarchical geometry orientation is applied. Each node
in the graph contains a rotation and translation with respect to its parent.
The translation is applied *first*. So if T is a translation operator, and
R is a rotation,

    x_final = Rx + p

In practice, both are performed using matrices.


Author: TJ Lane <tjlane@slac.stanford.edu>
June 11, 2015
"""

import re
import h5py
import warnings
import numpy as np
import scipy.ndimage.interpolation as interp

from psgeom import moveable
from psgeom import sensors
from psgeom import translate
from psgeom import basisgrid
from psgeom import metrology

_STRICT = False # global used for some testing purposes, ignore this


def arctan3(y, x):
    """
    Compute the inverse tangent. Like arctan2, but returns a value in [0,2pi].
    """
    theta = np.arctan2(y,x)
    if type(theta) == np.ndarray:
        theta[theta < 0.0] += 2 * np.pi
    else:
        if theta < 0.0: theta += 2 * np.pi
    return theta

     
class CompoundCamera(moveable.MoveableParent, moveable.MoveableObject):
    """
    The compound camera class contains its own local rotation and translation
    operations that provide a local frame for a set of children. The children
    can either be SensorElements or other CompoundCameras.
    """
    
    def __init__(self, type_name=None, id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):         
        """
        Create a CompoundCamera.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundCamera with the same name,
            if they are identical units. E.g., "QUAD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundCamera
            The parent frame, specified by an instance of CompoundCamera.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
        """
        
        self._type_name = type_name
        self._id        = id_num
        
        self.set_parent(parent)
        self._children = []
        
        self._rotation_angles = rotation_angles
        self._translation     = translation
        
        return
        
        
    @property
    def id_num(self):
        # todo this functionality is duplicated with self._id
        return self._id
        
    
    @property
    def num_pixels(self):
        return np.sum([ c.num_pixels for c in self._children ])
        
        
    @property
    def xyz(self):
        return np.array([ c.xyz for c in self._children ])


    def to_psana_file(self, filename, title='geometry'):
        """
        Write a geometry in psana format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Optional Parameters
        -------------------
        title : str
            Title of the geometry saved inside the file
        """
        translate.write_psana(self, filename, title)
        return
        
        
    @classmethod
    def from_psana_file(cls, filename):
        """
        Load a geometry in psana format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        root : CompoundCamera
            The CompoundCamera instance
            
        References
        ----------
        ..[1] https://confluence.slac.stanford.edu/display/PSDM/Detector+Geometry
        """
        ret = translate.load_psana(cls, filename)
        ret._sort_tree()
        return ret
            
    
class CompoundAreaCamera(CompoundCamera):
    """
    A specific kind of CompoundCamera, one with sensor elements that are
    planar rectangles. Most detectors should be CompoundAreaCameras.
    """
    
    def to_text_file(self, filename):
        """
        Write a geometry in raw text psf format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.
        """
        translate.write_psf_text(self, filename)
        return
        
    
    @classmethod
    def from_text_file(cls, filename):
        """
        Load a geometry in raw-text psf format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        root : detector.CompoundCamera
            The CompoundCamera instance
        """
        raise NotImplementedError()


    def to_hdf5(self, filename):
        """
        Save a geometry's xyz coordinates (self.xyz) in an HDF file.

        Parameters
        ----------
        filename : str
            The path of the file on disk.
        """

        f = h5py.File(filename, 'w')
        f['xyz'] = self.xyz
        f.close()

        return


    def from_hdf5(self, filename):
        raise NotImplementedError()
    

    def to_basisgrid(self):
        """
        Convert this object to a BasisGrid object, which represents the camera
        geometry as a set of vectors specifying the slow-scan and fast-scan
        edges of a set of panels
    
        Returns
        -------
        bg : basisgrid.BasisGrid
            The basisgrid object.
        """
    
        bg = basisgrid.BasisGrid()
    
        for sensor in self.leaves:
            if not isinstance(sensor, sensors.PixelArraySensor):
                raise TypeError('basisgrid representation is only compatible '
                                'with detectors that are entirely comprised of '
                                'PixelArrayElements')
            
            for g in sensor.psf:               
                bg.add_grid(*g)
    
        return bg
    

    @classmethod
    def from_basisgrid(cls, bg, element_type=sensors.Mtrx, 
                       type_name='root', strict=True):
        """
        Convert a BasisGrid object to a CompoundCamera.
    
        Parameters
        ----------
        bg : basisgrid.BasisGrid
            The basisgrid object to convert.
        
        element_type : sensors.PixelArraySensor
            The SensorElement type to populate the camera with.
        
        Returns
        -------
        cd : CompoundCamera
            The compound camera instance.
        """
    
        if not isinstance(bg, basisgrid.BasisGrid):
            raise TypeError('`bg` argument must be instance of BasisGrid,'
                            ' got: %s' % type(bg))
    
        cd = cls(type_name=type_name, id_num=0, parent=None)
    
    
        # loop over grids, possibly skipping some due to gaps...
        grid_index  = 0
        panel_index = 0
        while grid_index < bg.num_grids:
            
            p, s, f, shape = bg.get_grid(grid_index)
        
            pixel_shape = (np.linalg.norm(s),
                           np.linalg.norm(f))
            
            # FixedArraySensors know their shape/pixel_shape before hand and
            # can therefore convert basisgrids to gapped sensor geometries,
            # while the "else" clause attempts to be general and allow inference
            # of these values
            if issubclass(element_type, sensors.FixedArraySensor):
                pas = element_type(id_num=panel_index,
                                   parent=cd)
            else:
                pas = element_type(shape,
                                   pixel_shape,
                                   id_num=panel_index,
                                   parent=cd)
        
        
            # >>> determine panel rotation
            us = s / pixel_shape[0] # unit vector
            uf = f / pixel_shape[1] # unit vector
            n  = np.cross(uf, -us)  # tested for orthog. in next fxn
            
            ra = moveable._angles_from_rotated_frame(uf, -us, n)
            
            # >>> determine panel translation
            dims = pas.dimensions
            
            # center (-pixel_shape is due to the "half pixel" added by p)
            tr = p + (dims[0]-pixel_shape[0]) * us / 2.0 +\
                     (dims[1]-pixel_shape[1]) * uf / 2.0 
            
            pas._rotation_angles = ra
            pas._translation     = tr
            
            panel_index += 1
                
            # >>> deal with gaps
            # if gapped, we need to assemble many bgs into one sensor
            if pas.num_gaps > 0:
                    
                # check s/f vectors are all the same
                for ig in range(pas.num_gaps * 2): # each gap > 2 panels
                    
                    p2, s2, f2, shape2 = bg.get_grid( grid_index + ig )
                    
                    if not np.all(s2 == s) and np.all(f2 == f):
                        if strict:
                            raise ValueError('BasisGrid `bg` not compatible '
                                             'with `element_type`: %s, s/f '
                                             'vectors of gapped grids do not '
                                             'match' % str(element_type))
                        else:
                            print('WARNING s/f vectors of gapped grids do '
                                  'not match! (grid: %d & %d)' % (grid_index, ig))
            
                grid_index += ig + 1 # increment bg counter
                    
            # if no gaps, just goto next grid
            else: 
                grid_index  += 1
    
        return cd
    
    
    def to_crystfel_file(self, filename, coffset=None):
        """
        Write a geometry to disk in CrystFEL format. Note that some fields
        will be written but left blank -- these are fields you probably should
        fill in before performing any computations in CrystFEL, but are 
        information that we have no handle on (e.g. detector gain).
        When coffset is not given, coffset is set to detector distance and
        and clen is set to zero.

        Thanks to Rick Kirian & Tom White for assistance with this function.

        Parameters
        ----------
        filname : str
            The name of file to write. Will end in '.geom'

        coffset: float
            Detector home position to sample distance in metres
        """
        translate.write_generic_crystfel(self, filename, coffset=coffset)
        return
        
        
    @classmethod
    def from_crystfel_file(cls, filename, element_type=sensors.Mtrx):
        """
        Load a geometry in crystfel format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        camera : CompoundCamera
            The instance
        """
        return translate.load_crystfel(cls, filename, element_type=element_type)
        

# ---- specific detector implementations ---------------------------------------

class Cspad(CompoundAreaCamera):
    """
    This is for a 'full' size CSPAD. The need for a specific CSPAD object is
    rather unfortunate, but necessitated by assumptions made by other software
    we want to interact with.
    """

    @classmethod
    def from_basisgrid(cls, bg, element_type=sensors.Cspad2x1):
        
        # we want the element_type kwarg to maintain interface homogeneity
        # if element_type != sensors.Cspad2x1:
        #     raise RuntimeError('trying to create a CSPAD without Cspad2x1 elements!')
        
        cspad = super(Cspad, cls).from_basisgrid(bg, 
                                                 element_type=sensors.Cspad2x1, 
                                                 strict=False, 
                                                 type_name='CSPAD:V1')
        
        # enforce the CSPAD quad/asic heirarchy, which cannot be inferred from
        # the basisgrid alone
        for i, l in enumerate(cspad.leaves):
            if i % 8 == 0:
                quad = CompoundAreaCamera(type_name='QUAD:V1',
                                          id_num=i/8,
                                          parent=cspad)
            l.set_parent(quad)
                                
        return cspad
            

    def to_crystfel_file(self, filename, coffset=None, **kwargs):
        """
        Write a geometry to disk in CrystFEL format. Note that some fields
        will be written but left blank -- these are fields you probably should
        fill in before performing any computations in CrystFEL, but are 
        information that we have no handle on (e.g. detector gain).
        When coffset is not given, coffset is set to detector distance and
        and clen is set to zero.

        Thanks to Rick Kirian & Tom White for assistance with this function.

        Parameters
        ----------
        filname : str
            The name of file to write. Will end in '.geom'

        coffset: float
            Detector home position to sample distance in metres

        Optional Parameters
        -------------------
        maskfile : str
            Hdf5 filename of a mask used to indexing and integration by CrystFEL.
        """
        translate.write_cspad_crystfel(self, filename, coffset, 
                                       intensity_file_type='cheetah', 
                                       **kwargs)

        return
        
        
    def to_cheetah_file(self, filename):
        """
        Convert a CompoundCamera object to a Cheetah h5 pixel map.

        The CompoundCamera must be the a (4,8,185,388) shape representing a
        CSPAD -- cheetah can only understand CSPAD geometries. 

        Parameters
        ----------
        geometry : cspad.CSPad
            The detector geometry to write to disk
        filename : string
            The file name for the output pixel map
        """
        translate.write_cheetah(self, filename)
        return


    @classmethod
    def from_cheetah_file(cls, filename, element_type=sensors.Cspad2x1):
        """
        Load a geometry in cheetah format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        cspad : Cspad
            The Cspad instance
        """
        
        # we want the element_type kwarg to maintain interface homogeneity
        # if element_type != sensors.Cspad2x1:
        #     raise RuntimeError('trying to create a CSPAD without Cspad2x1 elements!')
            
        return translate.load_cheetah(cls, filename, element_type=sensors.Cspad2x1)


    @classmethod
    def from_metrology_file(cls, filename):
        """
        Load a geometry in metrology format.

        Note that the LCLS detector group often provides
        MS Excel files, but that this function expects a
        flat text, space delimited file of the form:

            # quad 0
            1 x1 y1 z1
            2 x2 y2 z2
            3 ...
            
            # quad 1
            1 x1 y1 z1
            2 x2 y2 z2
            3 ...

        Lines preceeded by a '#' will be ignored.

        It is recommened you simply generate this file by
        hand from whatever metrology information is provided,
        as historically there has not been a standard format
        as of the time of writing this (Sept 2018).

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        cspad : Cspad
            The Cspad instance

        """
        bg = metrology.load_to_basisgrid(filename)
        return cls.from_basisgrid(bg)



def load(filename, base=CompoundAreaCamera, infer_base=True):
    """
    Load a saved area camera from disk, attempting to interpert the 
    format from the file extension.

    Parameters
    ----------
    filename : str
        A path to the camera file on disk.

    base : camera.CompoundCamera
        The *class* (NOT instance!) that will be used to form
        the camera object. Must inheret camera.CompoundCamera.

    infer_base : bool
        If True, attempts to infer the base class (see `base`) 
        from the file header (#).

    Returns
    -------
    camera : camera.CompoundCamera
        The loaded camera object.
    """

    if infer_base:
        with open(filename, 'r') as f:
            text = f.read()

            # check for cspad
            to_find = re.compile("CSPAD|Cspad|CsPad|cspad")
            match_obj = to_find.search(text)
            if match_obj is not None:
                print(('Found `%s` in file, '
                      'interpreting geometry as CSPAD' % match_obj.group()))
                base = Cspad
    

    if not issubclass(base, CompoundCamera):
        raise TypeError('`base` of type %s is not an instance of CompoundCamera'
                        '' % type(base))


    if filename.endswith('.data'):
        camera = base.from_psana_file(filename)
    elif filename.endswith('.txt'):
        camera = base.from_text_file(filename)
    elif filename.endswith('.geom'):
        camera = base.from_crystfel_file(filename)
    elif filename.endswith('.h5'):
        camera = base.from_cheetah_file(filename)
    else:
        ext = filename.split('.')[-1]
        raise IOError('Could not understand extension: %s' % ext)


    return camera
