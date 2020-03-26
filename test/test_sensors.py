# ---- sensors.py --------------------------------------------------------------
#    + moveable.py

import numpy as np
import os
import unittest
import h5py

from psgeom import sensors
from psgeom import moveable
from psgeom import camera

class TestPixelArraySensor(object):
    
    def setup(self):
        
        self.shape = (185, 388)
        self.pixel_shape = (1.0, 1.0)
        
        self.rotation_angles = np.random.rand(3) * 360.0
        self.translation     = np.random.randn(3)
        
        self.PAS = sensors.PixelArraySensor(self.shape, self.pixel_shape, 
                                            type_name='Test', id_num=0, parent=None,
                                            rotation_angles=self.rotation_angles, 
                                            translation=self.translation)
    
    
    def test_local_transform(self):
        R = moveable._rotation_matrix_from_angles(*self.rotation_angles, 
                                         dummy_dimension=True)
        T = moveable._translation_matrix_from_vector(self.translation)
        L_ref = np.dot(T, R)
        L_ans = self.PAS.local_transform
        np.testing.assert_array_almost_equal(L_ref, L_ans)
    
    
    def test_global_transform(self):
        
        ra_p = np.random.rand(3) * 360.0
        t_p  = np.random.randn(3)
        
        Rp = moveable._rotation_matrix_from_angles(*ra_p, dummy_dimension=True)
        Tp = moveable._translation_matrix_from_vector(t_p)
        
        parent_obj = camera.CompoundCamera(type_name='daddy', 
                                               id_num=0, parent=None,
                                               rotation_angles=ra_p, 
                                               translation=t_p)
                                      
        self.PAS.set_parent(parent_obj)
        
        R = moveable._rotation_matrix_from_angles(*self.rotation_angles, 
                                         dummy_dimension=True)
        T = moveable._translation_matrix_from_vector(self.translation)
        
        G_ref = np.dot( np.dot( np.dot(Tp, Rp), T), R) # T_1 . R_1 . T_2 . R_2
        G_ans = self.PAS.global_transform
        
        np.testing.assert_array_almost_equal(G_ref, G_ans)
    
    
    def test_evaluate_transform(self):
    
        # for rotation
        x = np.array([[1.0, 0.0, 0.0]]) # vector pointing at x
        Rz = moveable._rotation_matrix_from_angles(90.0, 0.0, 0.0, 
                                                   dummy_dimension=True)                                         
        ref = np.array([[0.0, 1.0, 0.0]])
        ans = camera.CompoundCamera._evaluate_transform(Rz, x)
        np.testing.assert_array_almost_equal(ans, ref, err_msg='rotation')
                                         
        # for translation
        x = np.random.randint(0,5,size=(1,3))    
        y = np.random.randint(0,5,size=(1,3))    
        T = moveable._translation_matrix_from_vector(x)    
        assert np.all( camera.CompoundCamera._evaluate_transform(T, y) == x + y )
    
    
    def test_untransformed_xyz(self):
        uxyz = self.PAS.untransformed_xyz
        assert uxyz.shape[:-1] == self.shape, '%s %s' % (uxyz.shape, self.shape)
        # remember, slow is x convention... (changed 8-28-16 by TJL)
        assert np.max(uxyz[...,0]) == (self.shape[0]-1) * self.pixel_shape[0]
        assert np.max(uxyz[...,1]) == (self.shape[1]-1) * self.pixel_shape[1]
        
    
    def test_xyz(self):
        
        uxyz = self.PAS.untransformed_xyz
        
        buff = np.ones( list(uxyz.shape[:-1]) + [1], dtype=uxyz.dtype)
        uxyzd = np.concatenate([uxyz, buff], axis=-1)
        
        R = moveable._rotation_matrix_from_angles(*self.rotation_angles, 
                                                  dummy_dimension=True)
        T = moveable._translation_matrix_from_vector(self.translation)
        
        xyz_ans = np.dot(uxyzd, np.dot(T, R).T)
        np.testing.assert_array_almost_equal(self.PAS.xyz, xyz_ans[...,:3])
    

class TestSens2x1(TestPixelArraySensor):
    
    def test_2x1_central_gap(self):
        
        # regression test for the size of the big pixels
        
        s = sensors.Cspad2x1()
        
        small1 = s.xyz[0,193,0] - s.xyz[0,192,0]
        big    = s.xyz[0,194,0] - s.xyz[0,193,0]
        small2 = s.xyz[0,195,0] - s.xyz[0,194,0]
        
        np.testing.assert_array_almost_equal(small1, 109.92) # px size
        np.testing.assert_array_almost_equal(small2, 109.92)
        np.testing.assert_array_almost_equal(big,    439.68) # gap size     