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
in the graph contains a rotation and translation with respect to it's parent.
The translation is applied *first*. So if T is a translation operator, and
R is a rotation,

    x_final = Rx + p

In practice, both are performed using matrices.


Author: TJ Lane <tjlane@slac.stanford.edu>
June 11, 2015
"""

import re
import warnings
import numpy as np

from psgeom import moveable
from psgeom import sensors
from psgeom import translate
from psgeom import basisgrid

_STRICT = False # global used for some testing purposes, ignore this

     
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
        
        
    def add_child(self, child):
        """
        Add a child to the compound camera. This can either be a
        `SensorElement` or another `CompoundCamera`.
        
        Parameters
        ----------
        child : SensorElement or CompoundCamera
            The child object to add to this CompoundCamera node.
        """
        
        if not (isinstance(child, CompoundCamera) or \
                isinstance(child, sensors.SensorElement)):
                raise TypeError('`child` must be type: SensorElement or '
                                'CompoundCamera')
        
        for c in self.children:
            if c.name == child.name:
                if c is child:
                    raise NameError('Child object already registered with parent!')
                else:
                    raise NameError('Child with name %s already registered with'
                                    ' this parent (%s) -- please change the ID'
                                    ' number to give this object a unique name '
                                    'and re-register it as a child object' % \
                                    (child.name, self.name))
        
        self.children.append(child)
        child._parent = self
        
        return
        
        
    def draw_tree(self):
        """
        Sketch the camera tree, with this node as the root (higher levels in
        the heirarchy will not be shown)
        """

        print "--- " + str(self.name)
        
        def draw_child_tree(current, depth):
        
            for c in current.children:
                print depth * "    " + "|-- " + str(c.name)
                if hasattr(c, 'children'):
                    draw_child_tree(c, depth + 1)
                    
        draw_child_tree(self, 1)
        
        return
        
        
    def _sort_tree(self):
        """
        Order the tree by the id_num of each tree node.
        """
        
        self._children = sorted(self._children, key=lambda x : x.id_num)
        for c in self.children:
            if hasattr(c, '_sort_tree'):
                c._sort_tree()
        
        return
        
    
    @property
    def id_num(self):
        return self._id
    
        
    @property
    def num_children(self):
        return len(self._children)
    
        
    @property
    def children(self):
        return self._children
        
        
    @property
    def leaves(self):
        
        leaves = []
        
        def add_leaves(node):
            for c in node.children:
                if hasattr(c, 'children'):
                    add_leaves(c)
                else:
                    leaves.append(c)
                    
        add_leaves(self)
                    
        return leaves
        
    
    @property
    def num_pixels(self):
        return np.sum([ c.num_pixels for c in self._children ])
        
        
    @property
    def xyz(self):
        return np.array([ c.xyz for c in self._children ])


    def to_psana_file(self, filename, dist=1.0, title='geometry'):
        """
        Write a geometry in psana format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.
        
        dist : float
            Detector distance in metres.

        Optional Parameters
        -------------------
        title : str
            Title of the geometry saved inside the file
        """
        translate.write_psana(self, filename, dist, title)
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
                           
            p, s, f = sensor.psf 
            bg.add_grid(p, s, f, sensor.shape)
    
        return bg
    

    @classmethod
    def from_basisgrid(cls, bg, element_type=sensors.PixelArraySensor):
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
    
        cd = cls(type_name='root_frame', id_num=0, parent=None)
    
    
        for g in range(bg.num_grids):
        
            p, s, f, shape = bg.get_grid(g)
        
            pixel_shape = (np.linalg.norm(s),
                           np.linalg.norm(f))
        
            # to compute the rotation, find the 
            us = s / pixel_shape[0] # unit vector
            uf = f / pixel_shape[1] # unit vector
            n  = np.cross(us, uf)   # tested for orthog. in next fxn
        
            # remember: in the matrix convention (Mikhail uses), +x is slow
            # and +y is fast
            ra = moveable._angles_from_rotated_frame(us, uf, n)

            # translation is just p
            tr = p
        
            pas = element_type(shape, 
                               pixel_shape, 
                               type_name='grid_element_%dx%d' % shape, 
                               id_num=g, 
                               parent=cd,
                               rotation_angles=ra, 
                               translation=tr)

    
        return cd
    
    
    def to_crystfel_file(self, filename):
        """
        Write a geometry to disk in CrystFEL format. Note that some fields
        will be written but left blank -- these are fields you probably should
        fill in before performing any computations in CrystFEL, but are 
        information that we have no handle on (e.g. detector gain).

        Thanks to Rick Kirian & Tom White for assistance with this function.

        Parameters
        ----------
        filname : str
            The name of file to write. Will end in '.geom'
        """
        translate.write_generic_crystfel(self, filename)
        return
        
        
    @classmethod
    def from_crystfel_file(cls, filename):
        """
        Load a geometry in crystfel format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        cspad : Cspad
            The Cspad instance
        """
        return translate.load_crystfel(cls, filename)
        

# ---- specific detector implementations ---------------------------------------

class Cspad(CompoundAreaCamera):
    """
    This is for a 'full' size CSPAD. The need for a specific CSPAD object is
    rather unfortunate, but necessitated by assumptions made by other software
    we want to interact with.
    """
            
    
    @classmethod
    def from_basisgrid(cls, bg):
        """
        Convert a BasisGrid object to a Cspad. The BasisGrid must have 64 grids/
        panels, one for each CSPAD ASIC.
        
        Parameters
        ----------
        bg : basisgrid.BasisGrid
            The basisgrid object to convert.
            
        Returns
        -------
        cspad : Cspad
            The Cspad instance.
        """
        
        cspad = cls(type_name='CSPAD:V1')
        
        if not isinstance(bg, basisgrid.BasisGrid):
            raise TypeError('`bg` argument must be instance of BasisGrid,'
                            ' got: %s' % type(bg))
        
        
        # if the grids are asics, we can strip out every other one and then
        # treat them like 2x1s        
        if bg.num_pixels / bg.num_grids == 185 * 388: # two-by-one grids
            stride = 1
            
        elif bg.num_pixels / bg.num_grids == 185 * 194: # asic grids
            stride = 2
            
        else:
            raise RuntimeError('Could not tell if BasisGrid grid elements are '
                               'CSPAD 2x1s or ASICs. Pixels per element:'
                               '%d' % (bg.num_pixels / bg.num_grids,))
        
        for g in range(0, bg.num_grids, stride):
            
            asic_id = (g/stride) # index from 0 to 7
            
            # find the quad geometry
            quad_index = asic_id / 8
            
            # we just put the quads in a zero'd frame, no knowledge of absolute
            # orientations
            quad_rot = np.array([0.0, 0.0, 0.0])
            quad_trs = np.array([0.0, 0.0, 0.0])
                                 
                                 
            # add a new quad if necessary
            if asic_id % 8 == 0:
                quad = CompoundCamera(type_name='QUAD:V1',
                                        id_num=quad_index,
                                        parent=cspad,
                                        rotation_angles=quad_rot, 
                                        translation=quad_trs)
                                 
            p, s, f, shape = bg.get_grid(g)
            
            pixel_shape = (np.linalg.norm(s),
                           np.linalg.norm(f))
            
            
            # compute the rotation based on s/f vectors
            us = s / pixel_shape[0] # unit vector
            uf = f / pixel_shape[1] # unit vector
            
            # BIG WARNING: There is a minus sign on `us` below, which is needed
            # to account for the fact that the slow scan is swapped in the
            # sensors.Cspad2x1 class
            
            n  = np.cross(uf, -us)   # tested for orthog. in next fxn
            ra = moveable._angles_from_rotated_frame(uf, -us, n)
            
            
            # translation is center of 2x1, less the quad center
            # dont forget the big pixels!
            # dont forget p/s/f here is only an ASIC!
            # dont forget p points to the middle of a pixel!

            tr = p + 184.0/2.0 * 109.92 * us + (192.5 * 109.92 + 274.8) * uf

            
            # construct the 2x1
            pas = sensors.Cspad2x1(type_name='SENS2X1:V1', 
                                   id_num=int(asic_id % 8), 
                                   parent=quad,
                                   rotation_angles=ra, 
                                   translation=tr)
                        
                                   
            # if we skipped a grid (ASIC), we'll check to make sure that
            # the p-vector from that grid points to the correct pixel on
            # our newly oriented 2x1 -- allow 10 um error in x/y, 200 um in z
            
            if stride == 2:
                
                p_skipped, _, _, _ = bg.get_grid(g + 1)
                
                # be more lenient in z, since some programs are not general
                # enough to handle it
                
                if (np.linalg.norm(p_skipped[:2] - pas.xyz[0,194,:2]) > 10.0) or \
                   (np.abs(p_skipped[2] - pas.xyz[0,194,2]) > 200.0):
                    
                    print 'quad %d / 2x1 %d' % (quad_index, asic_id % 8)
                    print 'grid p-vector:   ', p_skipped
                    print 'pixel (0, 194):  ', pas.xyz[0,194,:]
                    print ''
                    
                    warnings.warn('The two ASICs making up the %d-th 2x1 on '
                                  'the %d-th quad (grids %d, %d) do not conform'
                                  ' to the geometric requirements of a 2x1 '
                                  'unit. Check your geometry! Do not ignore this'
                                  ' warning unless you were expecting it!!'
                                  '' % (asic_id % 8, quad_index, g, g+1))
                                  
                    if _STRICT:
                        raise RuntimeError('_STRICT set, no warnings allowed')
                               
                               
        return cspad
    
    
    def to_basisgrid(self):
        """
        Convet to a basisgrid where the individual ASICs are their own grids.
        
        Returns
        -------
        bg : basisgrid.BasisGrid
            The basisgrid object.
        """
        
        bg = basisgrid.BasisGrid()
        asic_shape = (185, 194)
        
        for sensor in self.leaves:
            if not isinstance(sensor, sensors.Cspad2x1):
                raise TypeError('basisgrid representation is only compatible '
                                'with detectors that are entirely comprised of '
                                'PixelArrayElements')
                               
            p, s, f = sensor.psf 
            
            # add the first ASIC of a 2x1...
            bg.add_grid(p, s, f, (185, 194))
            
            # then translate along the fast-scan dimension and add the second
            # DONT FORGET THE BIG PIXELS!!! (+3 pixels for gap)
            
            bg.add_grid(p + f * 197, s, f, (185, 194))
        
        return bg
    
        
    def to_crystfel_file(self, filename, coffset=0.0):
        """
        Write a geometry to disk in CrystFEL format. Note that some fields
        will be written but left blank -- these are fields you probably should
        fill in before performing any computations in CrystFEL, but are 
        information that we have no handle on (e.g. detector gain).

        Thanks to Rick Kirian & Tom White for assistance with this function.

        Parameters
        ----------
        filname : str
            The name of file to write. Will end in '.geom'

        coffset: float
            Detector home position to sample distance in metres
        """
        translate.write_cspad_crystfel(self, filename, coffset, intensity_file_type='cheetah')
        return
        
        
    @classmethod
    def from_crystfel_file(cls, filename):
        """
        Load a geometry in crystfel format.

        Parameters
        ----------
        filename : str
            The path of the file on disk.

        Returns
        -------
        cspad : Cspad
            The Cspad instance
        """
        return translate.load_crystfel(cls, filename)
        
        
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
    def from_cheetah_file(cls, filename):
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
        return translate.load_cheetah(cls, filename)
        


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
                print('Found `%s` in file, '
                      'interpreting geometry as CSPAD' % match_obj.group())
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
