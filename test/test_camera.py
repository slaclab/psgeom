# ---- camera.py -------------------------------------------------------------

import numpy as np
import os
import unittest
import h5py

from psgeom import translate
from psgeom import camera
from psgeom import sensors

PIXEL_TOLERANCE_um = 10.0

class TestCompoundCamera(object):
    
    def setup(self):
        self.geom = camera.CompoundCamera.from_psana_file('ref_files/cspad/refgeom_psana.data')
        
    
    def test_read_write(self):

        self.geom.to_psana_file('test.data')
        geom2 = camera.CompoundCamera.from_psana_file('test.data')
        
        assert self.geom.xyz.shape == geom2.xyz.shape, 'shape/element mismatch'            
        np.testing.assert_allclose(self.geom.xyz, geom2.xyz, rtol=1e-3)

        os.remove('test.data')
    
    def test_xyz_vs_old_implementation(self):
    
        # ---- get the geometry Mikhail-style
        #try:
        #    from PSCalib.GeometryAccess import GeometryAccess
        #    ga = GeometryAccess('ref_files/cspad/refgeom_psana.data')
        #    ref_xyz = ga.get_pixel_coords()
        #except:

        # if that don't work, load a pre-saved answer
        print('could not use GeometryAccess, loading saved xyz')
        ref_xyz = np.load('ref_files/cspad/GA_saved_1-end.npy')
        
        # some np-foo to move the 3-d x,y,z axis from first dim to last
        ref_xyz = np.rollaxis(np.array(ref_xyz), 0, 7) # send 0 --> 7
        ref_xyz = np.squeeze(ref_xyz)
    
        geom = camera.CompoundCamera.from_psana_file('ref_files/cspad/refgeom_psana.data')
        xyz_new = np.squeeze(geom.xyz)
    
        assert xyz_new.shape == ref_xyz.shape, 'shape mismatch'
    
        err = np.sum( np.abs(xyz_new - ref_xyz) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - ref_xyz) > 109.92 )
    
        assert err < 10.0, 'error greater than 10um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'
    
        #np.testing.assert_allclose(xyz_new, ref_xyz, atol=2 * 110.0)
        
        
    def test_leaves(self):
        assert len(self.geom.leaves) == 32, 'got: %d' % len(self.geom.leaves)
        
        
        
class TestCompoundAreaCamera:
    
    def setup(self):
        self.geom = camera.CompoundAreaCamera.from_psana_file('ref_files/cspad/refgeom_psana.data')
        self.klass = camera.CompoundAreaCamera


    def test_hdf5_file(self):
        self.geom.to_hdf5('ref_files/tmp_hdf.h5')

        f = h5py.File('ref_files/tmp_hdf.h5', 'r')
        xyz2 = np.array(f['/xyz'])
        f.close()
        np.testing.assert_allclose(self.geom.xyz, xyz2)
        os.remove('ref_files/tmp_hdf.h5')


    def test_rayonix_big(self):
        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/rayonix/big_rayonix.data')
        geom.to_crystfel_file('ref_files/tmp_rayonix.geom')
        os.remove('ref_files/tmp_rayonix.geom')


    def test_jungfrau_vs_geometry_access(self):

        ref_xyz = np.load('ref_files/jungfrau/jungfrau_saved.npy')

        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/jungfrau/jungfrau4m.data')
        xyz_new = np.squeeze(geom.xyz)

        assert xyz_new.shape == ref_xyz.shape, 'shape mismatch'
        print((xyz_new - ref_xyz)[0])

        err = np.sum( np.abs(xyz_new - ref_xyz) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - ref_xyz) > 75.0 )

        assert err < 10.0, 'error greater than 10um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'


    def test_bg_index_to_camera_index(self):
        # we have a cspad 2M: 4x4=16 leaves, 2 subpanels/leaf
    
        assert self.geom._bg_index_to_camera_index(0)  == (0,0)
        assert self.geom._bg_index_to_camera_index(1)  == (0,1)
        assert self.geom._bg_index_to_camera_index(2)  == (1,0)
        assert self.geom._bg_index_to_camera_index(3)  == (1,1)
        assert self.geom._bg_index_to_camera_index(31) == (15,1)

    
class TestCspad:
    
    def setup(self):

        self.geom = camera.CompoundAreaCamera.from_psana_file('ref_files/cspad/refgeom_psana.data')
        self.klass = camera.CompoundAreaCamera
        
    def test_to_basis_grid(self):

        bg = self.geom.to_basisgrid()
        xyz = np.squeeze(self.geom.xyz)

        for i in range(4):
            for j in range(8):

                bg_xyz_ij_1 = bg.grid_as_explicit(i*16 + j*2)
                bg_xyz_ij_2 = bg.grid_as_explicit(i*16 + j*2 + 1)
                
                cd_xyz_ij = xyz[i,j,:,:,:]
                
                np.testing.assert_allclose(bg_xyz_ij_1, cd_xyz_ij[:,:194,:], 
                                           atol=PIXEL_TOLERANCE_um)
                np.testing.assert_allclose(bg_xyz_ij_2, cd_xyz_ij[:,194:,:], 
                                           atol=PIXEL_TOLERANCE_um)
                                           
    
    def test_basisgrid_roundtrip(self):

        bg = self.geom.to_basisgrid()
        new = self.klass.from_basisgrid(bg, element_type=sensors.Cspad2x1)

        ref_xyz = np.squeeze(self.geom.xyz)
        new_xyz = new.xyz.reshape(ref_xyz.shape)

        assert self.geom.num_pixels == new.num_pixels
        
        np.testing.assert_allclose(ref_xyz, 
                                   new_xyz,
                                   atol=PIXEL_TOLERANCE_um)
        

    def test_hdf5_file(self):
        self.geom.to_hdf5('ref_files/tmp_hdf.h5')

        f = h5py.File('ref_files/tmp_hdf.h5', 'r')
        xyz2 = np.array(f['/xyz'])
        f.close()

        # we re-shape the xyz to match the way psana
        # presents CSPAD data
        assert xyz2.shape == self.geom.xyz.shape
        np.testing.assert_allclose(self.geom.xyz, xyz2)

        os.remove('ref_files/tmp_hdf.h5')
        
        
class TestJungfrau:
    
    def setup(self):
        
        self.geom = camera.CompoundAreaCamera.from_psana_file('ref_files/jungfrau/jungfrau4m.data')
        self.klass = camera.CompoundAreaCamera

    def test_bg_index_to_camera_index(self):
        # we have a JF 4M: 4x2=8 leaves, 8 subpanels/leaf

        for i in range(8):
            assert self.geom._bg_index_to_camera_index(i) == (0,i)
        assert self.geom._bg_index_to_camera_index(8)  == (1,0)
        assert self.geom._bg_index_to_camera_index(63) == (7,7)

        
    def test_to_basis_grid(self):

        bg = self.geom.to_basisgrid()
        xyz = np.squeeze(self.geom.xyz)

        for i in range(8):
            
            cd_xyz_ij = xyz[i,:,:]

            for j in range(2):
                for k in range(4):
                    
                    bg_idx = i*8 + j*4 + k
                    bg_xyz_ij = bg.grid_as_explicit(bg_idx)
                    
                    #print(bg_idx, i, j, k, bg_xyz_ij[0,0] - cd_xyz_ij[256*j:256*(j+1),256*k:256*(k+1)][0,0])
                    
                    np.testing.assert_allclose(bg_xyz_ij,
                                               cd_xyz_ij[256*j:256*(j+1),256*k:256*(k+1)],
                                               atol=PIXEL_TOLERANCE_um)
                                           
    
    def test_basisgrid_roundtrip(self):

        bg = self.geom.to_basisgrid()
        assert bg.num_grids == 8 * 8
        new = self.klass.from_basisgrid(bg, element_type=sensors.JungfrauSegment)

        ref_xyz = self.geom.xyz
        new_xyz = new.xyz
        
        assert self.geom.num_pixels == new.num_pixels
        assert ref_xyz.shape == new_xyz.shape
        
        np.testing.assert_allclose(ref_xyz, 
                                   new_xyz,
                                   atol=PIXEL_TOLERANCE_um)
                                   
                                   
    def test_basisgrid_roundtrip_random(self):
        
        geom = camera.CompoundAreaCamera()
        rand_trans = np.random.randn(3)
        rand_rot   = 90.0 * np.random.rand(3)
        
        pas = sensors.JungfrauSegment(parent=geom, id_num=0,
                                      rotation_angles=rand_rot,
                                      translation=rand_trans)
                                      
        bg = geom.to_basisgrid()
        assert bg.num_grids == 8
        new = camera.CompoundAreaCamera.from_basisgrid(bg, element_type=sensors.JungfrauSegment)
                                              
        ref_xyz = geom.xyz
        new_xyz = new.xyz
        
        assert ref_xyz.shape == new_xyz.shape
        
        np.testing.assert_allclose(ref_xyz, 
                                   new_xyz,
                                   atol=PIXEL_TOLERANCE_um)                              
        

        
class TestRayonix:

    # effectively tests Mtrx as well

    def test_v1(self):

        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/rayonix/rayonix_mtrxv1.data')
        ref_xyz = np.load('ref_files/rayonix/rayonix_saved_mtrxv1.npy')
    
        ref_xyz = np.rollaxis(np.array(ref_xyz), 0, 5) # send 0 --> end
        ref_xyz = np.squeeze(ref_xyz)
        xyz_new = np.squeeze(geom.xyz)
    
        assert xyz_new.shape == ref_xyz.shape, 'shape mismatch %s / %s' % (xyz_new.shape, ref_xyz.shape)
    
        err = np.sum( np.abs(xyz_new - ref_xyz) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - ref_xyz) > 89.0 )
    
        print('new', xyz_new)
        print('old', ref_xyz)
    
        assert err < 10.0, 'error greater than 10 um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px' 


    #def test_v2(self):
    #    pass

