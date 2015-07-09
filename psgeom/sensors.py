

"""
sensors.py
"""

import abc
import numpy as np

from psgeom import moveable


# ---- abstract sensor class  --------------------------------------------------

class SensorElement(moveable.MoveableObject):
            
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
        
    
# ---- specific sensor implementations  ----------------------------------------
    
class PixelArraySensor(SensorElement):

    def __init__(self, shape, pixel_shape, 
                 type_name='None', id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):
        """
        fill me in
        
        shape is (slow, fast)
        """
        
        self._type_name = type_name
        self._id        = id_num
        
        self.set_parent(parent)
        
        self._rotation_angles = rotation_angles
        self._translation     = translation
        
        self.shape       = tuple(shape)
        self._pixel_shape = np.array(pixel_shape)
        
        return
    
        
    @property
    def num_pixels(self):
        return np.product(self.shape)
    

    @property
    def untransformed_xyz(self):
                
        # convention that x is the quickly varying dimension (fast), y is slow
        # and z is perpendicular to the sensor in the untransformed view
        
        xy = np.mgrid[0.0:float(self.shape[1]),0.0:float(self.shape[0])].T
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]
        
        z = np.zeros([self.shape[0], self.shape[1], 1])
        
        xyz = np.concatenate([xy, z], axis=-1)
        
        return xyz
    
        
    @property
    def psf(self):
        """
        """
        
        xyz = self.xyz
        
        p = xyz[0,0,:]
        s = xyz[1,0,:] - p
        f = xyz[0,1,:] - p
            
        return p, s, f


# ---- specific sensor implementations  ------------------------------------------------

class Cspad2x1(PixelArraySensor):
    """
    
    https://confluence.slac.stanford.edu/display/PSDM/CSPAD+Geometry+and+Alignment
    """
    
    def __init__(self, **kwargs):
                 
        shape = (185, 388)
        pixel_shape = np.array([109.92, 109.92])
                 
        super(Cspad2x1, self).__init__(shape, pixel_shape, **kwargs)
                                       
        return
        
        
    @property
    def untransformed_xyz(self):

        # convention that x is the quickly varying dimension (fast), y is slow
        # and z is perpendicular to the sensor in the untransformed view

        xy = np.mgrid[0.0:float(self.shape[1]),0.0:float(self.shape[0])].T
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]

        # swap the slow-scan dimension to remain consistent with psana
        # TJL -- think -- does this violate the spec? it is documented at least.
        
        xy[:,:,:] = xy[::-1,:,:]
        
        # the CSPAD's central pixels are bigger than usual along the x dim
        # normal pixels are 109.92 x 109.92 um, the middle two columns are
        # 109.92 x 274.8 um. By translating the 2nd ASIC, we get most of the
        # pixels right, but the central columns will be a bit off
        
        xy[:,193:,0] += 2.0 * (274.8 - 109.92)        

        # and, finally, for some reason M measures rotations from the
        # center of the 2x1 but the corner of the quad. So we center the
        # sensor elements
        
        xy[:,:,0] -= np.mean(xy[:,:,0])
        xy[:,:,1] -= np.mean(xy[:,:,1])

        z = np.zeros([self.shape[0], self.shape[1], 1])

        xyz = np.concatenate([xy, z], axis=-1)

        return xyz
        
        
# ------------------------------------------------------------------------------
# define a "type map" that maps a list of known object identifier strings to
# the corresponding types

type_map = {'SENS2X1:V1' : Cspad2x1,
            'SENS2X1'    : Cspad2x1}






