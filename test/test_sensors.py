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
        
        # remember, slow is y convention... (changed Mar 27, 2020 by TJL)
        #           and sensor is centered
        assert np.max(uxyz[...,1]) == (self.shape[0]-1) * self.pixel_shape[0]/2.0
        assert np.max(uxyz[...,0]) == (self.shape[1]-1) * self.pixel_shape[1]/2.0
        
    
    def test_xyz(self):
        
        uxyz = self.PAS.untransformed_xyz
        
        buff = np.ones( list(uxyz.shape[:-1]) + [1], dtype=uxyz.dtype)
        uxyzd = np.concatenate([uxyz, buff], axis=-1)
        
        R = moveable._rotation_matrix_from_angles(*self.rotation_angles, 
                                                  dummy_dimension=True)
        T = moveable._translation_matrix_from_vector(self.translation)
        
        xyz_ans = np.dot(uxyzd, np.dot(T, R).T)
        np.testing.assert_array_almost_equal(self.PAS.xyz, xyz_ans[...,:3])
      

class TestGaps:

    def setup(self):
        
        self.shape = (128, 128)
        self.pixel_shape = (1.0, 1.0)
        self.pas = sensors.PixelArraySensor(self.shape, self.pixel_shape, 
                                            type_name='Test', id_num=0, 
                                            parent=None)
                                            

    def test_psf(self):
        
        r = self.pas.psf[0]
        assert np.all( r[0] == np.array([-63.5,  63.5 ,    0. ]) ), r[0]
        assert np.all( r[1] == np.array([ 0., -1.,  0.]) ), r[1]
        assert np.all( r[2] == np.array([1., 0., 0.]) ), r[2]
        assert r[3] == self.shape, r[3]
        
    
    def test_gap_access(self):
        
        self.pas.add_gap(2.0, 32, 'slow') # size, loc, axis
        self.pas.add_gap(2.0, 64, 'slow') # size, loc, axis
        self.pas.add_gap(2.0, 32, 'fast') # size, loc, axis
        self.pas.add_gap(2.0, 8,  'slow') # size, loc, axis
        self.pas.add_gap(2.0, 35, 'fast') # size, loc, axis
        
        sgs = self.pas._slow_gaps
        fgs = self.pas._fast_gaps
        
        assert len(sgs) == 3
        assert len(fgs) == 2
        
        cp = self.pas.shape[0]
        for g in sgs:
            assert g.location < cp
            cp = g.location
            
        cp = self.pas.shape[1]
        for g in fgs:
            assert g.location < cp
            cp = g.location
        
        
    def test_dimensions(self):
        
        self.pas.add_gap(2.5, 32, 'slow') # size, loc, axis
        
        d = self.pas.dimensions
        assert d[0] == 128.0 + 2.5
        assert d[1] == 128.0
        
        shp2 = (64, 64)
        ps2  = (0.5, 0.5)
        pas2 = sensors.PixelArraySensor(shp2, ps2, 
                                        type_name='Test', id_num=0, 
                                        parent=None)
        pas2.add_gap(2.5, 3, 'fast')
        pas2.add_gap(2.5, 5, 'fast')
        
        d2 = pas2.dimensions
        assert d2[0] == 32.0
        assert d2[1] == 34.5
        
        
    
    def test_slow_gap(self):
        
        self.pas.add_gap(2.0, 32, 'slow') # size, loc, axis
        
        r0 = self.pas.psf[0]
        r1 = self.pas.psf[1]
        
        # check shapes
        assert r0[3] == (32, 128)
        assert r1[3] == (128-32, 128)
        
        # original vectors
        assert np.all( r0[0] == np.array([-63.5,  65.0 ,    0. ]) ), r0[0]
        assert np.all( r0[1] == np.array([ 0., -1.,  0.]) ), r0[1]
        assert np.all( r0[2] == np.array([1., 0., 0.]) ), r0[2]
        
        # gapped vectors
        assert np.all( r1[0] == np.array([-63.5,  31.0,    0. ]) ), r1[0]
        assert np.all( r1[1] == np.array([ 0., -1.,  0.]) ), r1[1]
        assert np.all( r1[2] == np.array([1., 0., 0.]) ), r1[2]
        
        
    def test_fast_gap(self):
        
        self.pas.add_gap(2.0, 32, 'fast') # size, loc, axis
        
        r0 = self.pas.psf[0]
        r1 = self.pas.psf[1]
        
        # check shapes
        assert r0[3] == (128, 32)
        assert r1[3] == (128, 128-32)
        
        # original vectors
        assert np.all( r0[0] == np.array([-65.0,  63.5 ,    0. ]) ), r0[0]
        assert np.all( r0[1] == np.array([ 0., -1.,  0.]) ), r0[1]
        assert np.all( r0[2] == np.array([1., 0., 0.]) ), r0[2]
        
        # gapped vectors
        assert np.all( r1[0] == np.array([-31.0,  63.5,    0. ]) ), r1[0]
        assert np.all( r1[1] == np.array([ 0., -1.,  0.]) ), r1[1]
        assert np.all( r1[2] == np.array([1., 0., 0.]) ), r1[2]
        
        
    def test_pdf_many_gaps(self):
        
        self.pas.add_gap(2.0, 32, 'fast') # size, loc, axis
        self.pas.add_gap(2.0, 32, 'slow') # size, loc, axis
        
        shapes = [ r[3] for r in self.pas.psf ]
    
        assert shapes[0] == (32,32)
        assert shapes[1] == (32,96)
        assert shapes[2] == (96,32)
        assert shapes[3] == (96,96)
        
        
    def test_trans_bg_to_sensor(self):
        
        data = np.random.randn(*self.shape)
        assert np.all(data == self.pas.trans_bg_to_sensor(data)) # no gaps
        
        self.pas.add_gap(1.0, 32, 'fast') # size, loc, axis
        self.pas.add_gap(2.0, 32, 'slow') # size, loc, axis
        bg_data = [ data[:32,:32], data[:32,32:], data[32:,:32], data[32:,32:] ]
        
        assert np.all(data == self.pas.trans_bg_to_sensor(bg_data))
        
        
    def test_trans_sensor_to_bg(self):
        
        data = np.random.randn(*self.shape)
        assert np.all(data == self.pas.trans_sensor_to_bg(data)) # no gaps
        
        self.pas.add_gap(2.0, 32, 'fast') # size, loc, axis
        self.pas.add_gap(1.5, 32, 'slow') # size, loc, axis
        
        td = self.pas.trans_sensor_to_bg(data)
        assert len(td) == 4
        
        assert td[0].shape == (32,32) # topleft
        assert np.all(td[0] == data[:32,:32])
        
        assert td[1].shape == (32,96) # topright
        assert np.all(td[1] == data[:32,32:])
        
        assert td[2].shape == (96,32) # bottomleft
        assert np.all(td[2] == data[32:,:32])
        
        assert td[3].shape == (96,96) # bottomright
        assert np.all(td[3] == data[32:,32:])
        
        # test multi-panel pass
        data2 = np.random.randn(*(6,) + self.shape)
        td2 = self.pas.trans_sensor_to_bg(data2)
        assert len(td2) == 6 * 4
        

    def test_trans_consistency(self):
        
        data = np.random.randn(*self.shape)
        
        self.pas.add_gap(2.0, 14,  'fast') # size, loc, axis
        self.pas.add_gap(1.5, 103, 'slow') # size, loc, axis
        
        s_data = self.pas.trans_sensor_to_bg(data)
        r_data = self.pas.trans_bg_to_sensor(s_data)
        assert np.all(data == r_data)
        

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