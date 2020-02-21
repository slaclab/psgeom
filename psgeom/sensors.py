

"""
sensors.py
"""

import abc
import numpy as np

from psgeom import moveable


# ---- abstract sensor class  --------------------------------------------------

class SensorElement(moveable.MoveableObject):
    """
    Abstract base class specifying a SensorElement. These elements are the
    actual detecting units of the camera (e.g. a two-by-one for a CSPAD). Many
    such elements usually form a camera.
    
    These objects specify the location of pixels within a sensing element and
    a few other basic facts about them.
    """
            
    @abc.abstractproperty
    def num_pixels(self):
        return

        
    @property
    def pixel_shape(self):
        return self._pixel_shape

        
    @abc.abstractproperty
    def untransformed_xyz(self):
        pass

    
    @property    
    def xyz(self):
        uxyz = self.untransformed_xyz
        T = self.global_transform
        return self._evaluate_transform(T, uxyz)
        
        
    @property
    def id_num(self):
        return self._id


    @classmethod
    def from_type(cls, type_name,
                  id_num=None,
                  parent=None,
                  rotation_angles=None,
                  translation=None):
        raise NotImplementedError('from_type method not implemented')
        return        
        
    
# ---- specific sensor implementations  ----------------------------------------
    
class PixelArraySensor(SensorElement):
    """
    The PixelArraySensor is an implementation of a rectangular array sensor.
    Likely most cameras can be generated from instances of this element.
    """

    def __init__(self, shape, pixel_shape, 
                 type_name='None', id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):
        """
        Create a PixelArraySensor.
        
        Parameters
        ----------
        shape : tuple of ints
            The shape of the rectangular array of pixels. The ordering of the
            shape is (slow, fast), where fast indicates the direction that is
            most rapidly scanned across when mapping intensities stored in
            a linear memory array onto the camera.
            
        pixel_shape : tuple of floats
            The size of the rectangular pixels that make up the detector, also
            in (slow, fast) directions. Note that units in psgeom are arbitrary,
            but that you need to be consistent!
        
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "QUAD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : PixelArraySensor
            The sensor element.
        """
        
        self._type_name = type_name
        self._id        = id_num
        
        self.set_parent(parent)
        
        self._rotation_angles = rotation_angles
        self._translation     = translation
        
        self.shape        = tuple(shape)
        self._pixel_shape = np.array(pixel_shape)
        
        return
    
        
    @property
    def num_pixels(self):
        return np.product(self.shape)
    

    @property
    def untransformed_xyz(self):
        """
        Return the xyz coordinates of the element in the reference frame, that
        is before any translation/rotation operations have been applied.
        """
                
        # convention that x/row/first-index is the SLOW varying dimension
        # y/column/second-index is FAST and z is perpendicular to the sensor 
        # completing a right handed coordinate system in the untransformed view
        
        xy = np.mgrid[0.0:float(self.shape[0]),0.0:float(self.shape[1])]
        xy = np.rollaxis(xy, 0, start=3)
        
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]
        
        # add the z dimension (just flat)
        z = np.zeros([self.shape[0], self.shape[1], 1])
        xyz = np.concatenate([xy, z], axis=-1)
        
        return xyz   
        
    @property
    def psf(self):
        """
        Return basis grids for this object.
        
        Returns
        -------
        p : np.ndarray
            A 3-vector pointing from the interaction site to the first pixel
            read out from memory for this element.
            
        s : np.ndarray
            A 3-vector pointing along the slow scan direction. The size of the
            vector is the size of the pixel in this direction.
            
        f : np.ndarray
            A 3-vector pointing along the slow scan direction. The size of the
            vector is the size of the pixel in this direction.
        """
        
        xyz = self.xyz
        
        p = xyz[0,0,:]
        s = xyz[1,0,:] - p
        f = xyz[0,1,:] - p
            
        return p, s, f


    @classmethod
    def from_type(cls, type_name,
                  id_num=0, parent=None,
                  rotation_angles=np.array([0.0, 0.0, 0.0]),
                  translation=np.array([0.0, 0.0, 0.0])):
        """
        Factory function for automatically identifying
        the sensor based on the `type_name` alone.            
        """
        # see code in translate.py
        return cls(type_name=type_name, id_num=id_num, parent=parent,
                   rotation_angles=rotation_angles,
                   translation=translation)


# ---- specific sensor implementations  ------------------------------------------------

class Cspad2x1(PixelArraySensor):
    """
    A specific PixelArraySensor representing a CSPAD 2x1.
    """
    
    def __init__(self, **kwargs):
        """
        Create a Cspad2x1.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "QUAD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : Cspad2x1
            The sensor element.
        """
                 
        shape = (185, 388)
        pixel_shape = np.array([109.92, 109.92])
                 
        super(Cspad2x1, self).__init__(shape, pixel_shape, **kwargs)
                                       
        return

    
    @property
    def untransformed_xyz(self):
        """
        Return the xyz coordinates of the element in the reference frame, that
        is before any translation/rotation operations have been applied.
        """

        # of course, the CSPAD layout is different
        xy = np.mgrid[0.0:float(self.shape[1]),0.0:float(self.shape[0])].T
        xy[:,:,:] = xy[::-1,:,:]
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]
        
        # add the z dimension (just flat)
        z = np.zeros([self.shape[0], self.shape[1], 1])
        xyz = np.concatenate([xy, z], axis=-1)
        
        # the CSPAD's central pixels are bigger than usual along the x dim
        # normal pixels are 109.92 x 109.92 um, the middle two columns are
        # 109.92 x 274.8 um. By translating the 2nd ASIC, we get most of the
        # pixels right, but the central columns will be a bit off
        
        # this is equivalent to a 3-pixel shift
        # note that 2 * (274.80 - 109.92) = 329.76
        # gap is between pixel indices 193 & 194
        
        xyz[:,194:,0] += 2.0 * (274.8 - 109.92)        

        # and, finally, for some reason M measures rotations from the
        # center of the 2x1 but the corner of the quad. So we center the
        # sensor elements
        
        xyz[:,:,0] -= np.mean(xyz[:,:,0])
        xyz[:,:,1] -= np.mean(xyz[:,:,1])

        return xyz
        

class Mtrx(PixelArraySensor):
    """
    A specific PixelArraySensor representing a Rayonix sensor element.
    """
    
    def __init__(self, shape, pixel_shape, id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):
        """
        Create a Mtrx.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "RYONIX:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : Mtrx
            The sensor element.
        """
                 
        if shape is None or pixel_shape is None:
            raise RuntimeError()

        # TJL 4/9/18
        # I am not sure why these lines are necessary
        # but they seem to be to get these attributes set
        # I would have expected the super init method below to take care of it...
        self.shape = shape
        self._pixel_shape = pixel_shape

        super(Mtrx, self).__init__(shape, pixel_shape, 
                 type_name='shouldbeoverwritten', 
                 id_num=id_num, parent=parent,
                 rotation_angles=rotation_angles, 
                 translation=translation)
                                       
        return


    @classmethod
    def from_type(cls, type_name,
                  id_num=0, parent=None,
                  rotation_angles=np.array([0.0, 0.0, 0.0]),
                  translation=np.array([0.0, 0.0, 0.0])):
        """
        Factory function for automatically identifying
        the sensor based on the `type_name` alone.            
        """

        if 'MTRX' not in type_name:
            raise ValueError('`type_name` (%s) does not contain "MTRX"'
                             'cannot generate Mtrx object from type_name'
                             ' alone')

        s0, s1, ps0, ps1 = type_name.split(':')[1:]
        shape = (int(s0), int(s1))
        pixel_shape = (float(ps0), float(ps1))
        #print "###: ", shape, pixel_shape

        return cls(shape, pixel_shape,
                   id_num=id_num, parent=parent,
                   rotation_angles=rotation_angles,
                   translation=translation)

    @property
    def type_name(self):
        return 'MTRX:%d:%d:%d:%d' %(self.shape[0], self.shape[1], 
                                    round(self._pixel_shape[0]),
                                    round(self._pixel_shape[1]))

        
class PnccdQuad(PixelArraySensor):
    """
    A specific PixelArraySensor representing a pnCCD quad.
    """
    
    def __init__(self, **kwargs):
        """
        Create a PnccdQuad.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "PNCCD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : PnccdQuad
            The sensor element.
        """
                 
        shape = (512, 512)
        pixel_shape = np.array([75.0, 75.0]) # microns
                 
        super(PnccdQuad, self).__init__(shape, pixel_shape, **kwargs)
                                       
        return


class JungfrauSegment(PixelArraySensor):
    """
    A specific PixelArraySensor representing a 1M JUNGFRAU segment.
    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "JUNGFRAU:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : JungfrauSegment
            The sensor element.
        """

        shape = (512, 1024)
        pixel_shape = np.array([75.0, 75.0]) # microns

        super(JungfrauSegment, self).__init__(shape, pixel_shape, **kwargs)

        return


    @property
    def untransformed_xyz(self):
        """
        Return the xyz coordinates of the element in the reference frame, that
        is before any translation/rotation operations have been applied.
        """

        # JUNGFRAU 1M SENSOR LOOKS LIKE THIS
        #
        # all lines are a 2 pixel gap
        #
        #        fast -> (axis 1)
        #     -------------------------------------------------
        #  s  |           |           |           |           |
        #  l  | 256 x 256 | 256 x 256 | 256 x 256 | 256 x 256 |
        #  o  |           |           |           |           |
        #  w  -------------------------------------------------
        #  |  |           |           |           |           |
        #  0  | 256 x 256 | 256 x 256 | 256 x 256 | 256 x 256 |
        #     |           |           |           |           |
        #     -------------------------------------------------
        

        xy = np.mgrid[0.0:float(self.shape[1]),0.0:float(self.shape[0])].T
        xy[:,:,:] = xy[::-1,:,:]
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]

        # add the z dimension (just flat)
        z = np.zeros([self.shape[0], self.shape[1], 1])
        xyz = np.concatenate([xy, z], axis=-1)

        # add the two pixel gaps
        gap_size = 2.0

        xyz[:,256:,0] += gap_size * self.pixel_shape[1]
        xyz[:,512:,0] += gap_size * self.pixel_shape[1]
        xyz[:,768:,0] += gap_size * self.pixel_shape[1]
        xyz[256:,:,1] -= gap_size * self.pixel_shape[0] # minus due to flip convention

        # and, finally, for some reason M measures rotations from the
        # center of the 2x1 but the corner of the quad. So we center the
        # sensor elements
        xyz[:,:,0] -= np.mean(xyz[:,:,0])
        xyz[:,:,1] -= np.mean(xyz[:,:,1])

        return xyz

