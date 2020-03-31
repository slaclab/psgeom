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

        print("--- " + str(self.name))
        
        def draw_child_tree(current, depth):
        
            for c in current.children:
                print(depth * "    " + "|-- " + str(c.name))
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
    def from_basisgrid(cls, bg, element_type=sensors.Mtrx):
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
            print('***', grid_index, p, tr)
            
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
                            raise ValueError('BasisGrid `bg` not compatible '
                                             'with `element_type`: %s, s/f '
                                             'vectors of gapped grids do not '
                                             'match' % str(element_type))
                
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
    def from_crystfel_file(cls, filename):
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
        return translate.load_crystfel(cls, filename)
        

# ---- specific detector implementations ---------------------------------------

class Cspad(CompoundAreaCamera):
    """
    This is for a 'full' size CSPAD. The need for a specific CSPAD object is
    rather unfortunate, but necessitated by assumptions made by other software
    we want to interact with.
    """


    def assemble_image(self, raw_image):
        """
        Visualize CSPAD data, assembled using the geometry.

        The final image generated is a 2D approximation to the
        image as seen looking DOWNSTREAM (from the source).

        Parameters
        ----------
        raw_image : np.ndarray
            A shape (32,185,388) array that contains the image
            to be visualized.
        """
        
        # set up the raw image and the assembled template
        if not raw_image.shape == (32,185,388):
            raise ValueError('`raw_image` must have shape (32,185,388), got '
                             '%s' % str(raw_image.shape))
        
        # for some reason, bool types don't work. Make them ints
        if raw_image.dtype == np.bool:
            raw_image = raw_image.astype(np.int32)
        
        bounds = 2000 # JAS: total image range is 2000, ensures beam center is at (1000,1000)
        assembled_image = np.zeros((bounds, bounds), dtype=raw_image.dtype)
    
        bg = self.to_basisgrid()

        # iterate over quads
        pixel_size = 109.920
        for quad_index in range(4):
            for two_by_one in range(8):

                asic_idx = quad_index * 16 + two_by_one * 2 # add one for 2nd asic
                
                # assemble the 2x1 -- insert a 3 px gap
                gap = np.zeros( (185,3), dtype=raw_image.dtype )
                two_by_one_img = np.hstack( (raw_image[quad_index*8+two_by_one,:,:194], gap, 
                                             raw_image[quad_index*8+two_by_one,:,194:]) )
                
                # flip x data to conform w/CXI convention
                #two_by_one_img = two_by_one_img[::-1,:]
                
                # note that which dim is x changes w/two_by_one and quad_index
                # here the rotation is off between dtc/cspad by 180 in some quads
                # JAS: updated rotation to asic_rot - 180 instead of -asic_rot 
                #      to get proper rotation of asics in assembled image
                p, s, f, shape = bg.get_grid(asic_idx)
                theta = arctan3(f[1], f[0]) * (360. / (np.pi * 2.0))

                two_by_one_img = interp.rotate(two_by_one_img,
                                               theta - 180,
                                               output=two_by_one_img.dtype,
                                               reshape=True)
                
                # find the center of the 2x1 in space
                corners0 = bg.get_grid_corners(asic_idx)
                corners1 = bg.get_grid_corners(asic_idx + 1)
                
                # un-swap x-axis and re-swap below -- necessary b/c now we
                # have data in two_by_one_img that needs swap
                corners0[:,0] = -corners0[:,0]
                corners1[:,0] = -corners1[:,0]
                
                center = ( np.concatenate([corners0[:,0], corners1[:,0]]).mean(),
                           np.concatenate([corners0[:,1], corners1[:,1]]).mean() )

                # find the bottom left corner (note x is cols, so swap inds)
                c = (center[0] / pixel_size - two_by_one_img.shape[1] / 2.,
                     center[1] / pixel_size - two_by_one_img.shape[0] / 2.,)
                
                # the assembled image center will be at 1000, 1000 by convention
                cs = int(round(c[0])) + 1000
                rs = int(round(c[1])) + 1000

                if (rs < 0) or (rs+two_by_one_img.shape[0] > bounds):
                    raise ValueError('rs: out of bounds in rows. CSPAD geometry '
                                     'extends beyond 2000 x 2000 grid it is '
                                     'assembled on. It is likely that your CSPAD '
                                     'geometry is wacky in some respect -- use '
                                     '`sketch` method to check.')
                if (cs < 0) or (cs+two_by_one_img.shape[1] > bounds):
                    raise ValueError('cs: out of bounds in cols. CSPAD geometry '
                                     'extends beyond 2000 x 2000 grid it is '
                                     'assembled on. It is likely that your CSPAD '
                                     'geometry is wacky in some respect -- use '
                                     '`sketch` method to check.')
                
                assembled_image[rs:rs+two_by_one_img.shape[0],
                                cs:cs+two_by_one_img.shape[1]] += two_by_one_img
        
        # swap x-axis to conform to CXI convention
        #assembled_image = assembled_image[:,::-1]
        
        return assembled_image


    def sketch(self, mpl_axes=None, quad_colors = ['k', 'g', 'purple', 'b']):
        """
        Draw a rough sketch of the layout of the CSPAD
        """

        pixel_positions = np.squeeze(self.xyz)
        print(pixel_positions.shape)
        
        if not mpl_axes:
            from matplotlib import pyplot as plt
            import matplotlib.patches as plt_patches
            plt.figure()
            ax = plt.subplot(111)
        else:
            ax = mpl_axes

        for i in range(4):
            for j in range(8):
                x = pixel_positions[i,j,:,:,0]
                y = pixel_positions[i,j,:,:,1]
                corners = np.zeros((5,2))

                corners[0,:] = np.array([ x[0,0],   y[0,0] ])     # bottom left
                corners[1,:] = np.array([ x[0,-1],  y[0,-1] ])    # bottom right
                corners[3,:] = np.array([ x[-1,0],  y[-1,0] ])    # top left
                corners[2,:] = np.array([ x[-1,-1], y[-1,-1] ])   # top right
                corners[4,:] = np.array([ x[0,0],   y[0,0] ])     # make rectangle

                ax.plot(corners[:,0], corners[:,1], lw=2, color=quad_colors[i])
                ax.scatter(x[0,0], y[0,0])
                
        beam_center = plt_patches.Circle((0, 0), 2, fill=True, lw=1, color='orange')
        ax.add_patch(beam_center)
                
        # mirror x axis for CXI convention
        if not ax.xaxis_inverted():
            ax.invert_xaxis()

        if mpl_axes:
            return ax
        else:
            plt.show()
            return


    def imshow_cspad(self, image, vmin=0, vmax=None, mpl_axes=None):
        """
        Show an assembled image (e.g. from CSPad(raw_image) ) as it would be seen
        when viewed from upstream at CXI. CXI convention is that the plus-x direction
        is towards the hutch door, plus-y is upwards, and plus-z is the direction
        of the beam.
        
        Parameters
        ----------
        image : np.ndarray
            A two-dimensional assembled image
        
        Returns
        -------
        im : axes.imshow
            The imshow instance.
        """

        if image.shape == (32, 185, 388):
            img = self.assemble_image(image)
        else:
            img = image

        # to be consistent with CXI convention, we want +x going left, and +y up
        if not mpl_axes:
            from matplotlib import pyplot as plt
            plt.figure()
            ax = plt.subplot(111)
        else:
            ax = mpl_axes

        im = ax.imshow( img, origin='lower', vmin=vmin, vmax=vmax,
                        interpolation='nearest' )

        if not ax.xaxis_inverted():
            ax.invert_xaxis()

        return im

    
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
                    
                    print('quad %d / 2x1 %d' % (quad_index, asic_id % 8))
                    print('grid p-vector:   ', p_skipped)
                    print('pixel (0, 194):  ', pas.xyz[0,194,:])
                    print('')
                    
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
                               
            p, s, f, shp = sensor.psf[0] # TODO
            
            # add the first ASIC of a 2x1...
            bg.add_grid(p, s, f, (185, 194))
            
            # then translate along the fast-scan dimension and add the second
            # DONT FORGET THE BIG PIXELS!!! (+3 pixels for gap)
            
            bg.add_grid(p + f * 197, s, f, (185, 194))
        
        return bg


    def to_hdf5(self, filename):
        """
        Save a geometry's xyz coordinates (self.xyz) in an HDF file.

        Parameters
        ----------
        filename : str
            The path of the file on disk.
        """

        f = h5py.File(filename, 'w')
        f['xyz'] = np.vstack(np.squeeze(self.xyz))
        f.close()

        return
    

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
        translate.write_cspad_crystfel(self, filename, coffset, intensity_file_type='cheetah', **kwargs)

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
