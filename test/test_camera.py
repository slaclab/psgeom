# ---- camera.py -------------------------------------------------------------

import numpy as np
import os
import unittest
import h5py

from psgeom import translate
from psgeom import camera

PIXEL_TOLERANCE_um = 10.0

class TestCompoundCamera(object):
    
    def setup(self):
        self.geom = camera.CompoundCamera.from_psana_file('ref_files/refgeom_psana.data')
        
    
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
        #    ga = GeometryAccess('ref_files/refgeom_psana.data')
        #    xyz_old = ga.get_pixel_coords()
        #except:

        # if that don't work, load a pre-saved answer
        print('could not use GeometryAccess, loading saved xyz')
        xyz_old = np.load('ref_files/GA_saved_1-end.npy')
        
        # some np-foo to move the 3-d x,y,z axis from first dim to last
        xyz_old = np.rollaxis(np.array(xyz_old), 0, 7) # send 0 --> 7
        xyz_old = np.squeeze(xyz_old)
    
        geom = camera.CompoundCamera.from_psana_file('ref_files/refgeom_psana.data')
        xyz_new = np.squeeze(geom.xyz)
    
        assert xyz_new.shape == xyz_old.shape, 'shape mismatch'
    
        err = np.sum( np.abs(xyz_new - xyz_old) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - xyz_old) > 109.92 )
    
        assert err < 10.0, 'error greater than 10um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'
    
        #np.testing.assert_allclose(xyz_new, xyz_old, atol=2 * 110.0)
        
        
    def test_leaves(self):
        assert len(self.geom.leaves) == 32, 'got: %d' % len(self.geom.leaves)
        
        
        
class TestCompoundAreaCamera(TestCompoundCamera):
    
    def setup(self):
        self.geom = camera.CompoundAreaCamera.from_psana_file('ref_files/refgeom_psana.data')
        self.klass = camera.CompoundAreaCamera


    def test_hdf5_file(self):
        self.geom.to_hdf5('ref_files/tmp_hdf.h5')

        f = h5py.File('ref_files/tmp_hdf.h5', 'r')
        xyz2 = np.array(f['/xyz'])
        f.close()
        np.testing.assert_allclose(self.geom.xyz, xyz2)
        os.remove('ref_files/tmp_hdf.h5')

        
    def test_rayonix_vs_geometry_access(self):
    
        # ---- get the geometry Mikhail-style
        #try:
        #    from PSCalib.GeometryAccess import GeometryAccess
        #    ga = GeometryAccess('ref_files/rayonix.data')
        #    xyz_old = ga.get_pixel_coords()
        #except:

        # if that don't work, load a pre-saved answer
        print('could not use GeometryAccess, loading saved xyz')
        xyz_old = np.load('ref_files/rayonix_saved.npy')
    
        xyz_old = np.rollaxis(np.array(xyz_old), 0, 5) # send 0 --> end
        xyz_old = np.squeeze(xyz_old)
    
        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/rayonix.data')
        xyz_new = np.squeeze(geom.xyz)
    
        assert xyz_new.shape == xyz_old.shape, 'shape mismatch %s / %s' % (xyz_new.shape, xyz_old.shape)
    
        err = np.sum( np.abs(xyz_new - xyz_old) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - xyz_old) > 89.0 )
    
        print('new', xyz_new)
        print('old', xyz_old)
    
        assert err < 10.0, 'error greater than 10 um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'
        
        
    def test_rayonix_crystfel(self):
        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/rayonix.data')
        geom.to_crystfel_file('ref_files/tmp_rayonix.geom')
        
        # compare to reference created by make_pixelmap
        f = h5py.File('ref_files/rayonix.h5')
        x = np.array(f['/x']) * 1e6
        y = np.array(f['/y']) * 1e6
        #z = np.array(f['/z']) # z currently is not populated by make_pixelmap
        f.close()
        
        # NOTE: x is flipped in the crystFEL/cheetah convention
        # as there, z points towards the source of x-rays (not along travel)
        np.testing.assert_allclose( np.squeeze(geom.xyz[...,0]), -x, atol=10 )
        np.testing.assert_allclose( np.squeeze(geom.xyz[...,1]),  y, atol=10 )
        
        geom2 = camera.CompoundAreaCamera.from_crystfel_file('ref_files/tmp_rayonix.geom')
        np.testing.assert_allclose(geom.xyz, geom2.xyz, atol=0.001)
        
        os.remove('ref_files/tmp_rayonix.geom')


    def test_rayonix_big(self):
        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/big_rayonix.data')
        geom.to_crystfel_file('ref_files/tmp_rayonix.geom')
        os.remove('ref_files/tmp_rayonix.geom')
        
    
    def test_pnccd_vs_geometry_access(self):

        # ---- get the geometry Mikhail-style
        #try:
        #    from PSCalib.GeometryAccess import GeometryAccess
        #    ga = GeometryAccess('ref_files/pnccd.data')
        #    xyz_old = ga.get_pixel_coords()
        #except:

        # if that don't work, load a pre-saved answer
        print('could not use GeometryAccess, loading saved xyz')
        xyz_old = np.load('ref_files/pnccd_saved.npy')

        xyz_old = np.rollaxis(np.array(xyz_old), 0, 6) # send 0 --> end
        xyz_old = np.squeeze(xyz_old)

        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/pnccd.data')
        xyz_new = np.squeeze(geom.xyz)

        assert xyz_new.shape == xyz_old.shape, 'shape mismatch'

        err = np.sum( np.abs(xyz_new - xyz_old) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - xyz_old) > 75.0 )

        assert err < 10.0, 'error greater than 10um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'


    def test_jungfrau_vs_geometry_access(self):

        xyz_old = np.load('ref_files/jungfrau_saved.npy')

        geom = camera.CompoundAreaCamera.from_psana_file('ref_files/refgeom_jungfrau4m.data')
        xyz_new = np.squeeze(geom.xyz)

        assert xyz_new.shape == xyz_old.shape, 'shape mismatch'
        print((xyz_new - xyz_old)[0])

        err = np.sum( np.abs(xyz_new - xyz_old) ) / float(np.product(xyz_new.shape))
        print('Mean Absolute Error: %f um / px' % err)
        num_more_than_1px_err = np.sum( np.abs(xyz_new - xyz_old) > 75.0 )

        assert err < 10.0, 'error greater than 10um avg per px (%f)' % err
        assert num_more_than_1px_err < 7500, '>7500 pix w err > 1 px'
    
    
class TestCspad(TestCompoundCamera):
    
    def setup(self):
        self.geom = camera.Cspad.from_psana_file('ref_files/refgeom_psana.data')
        self.klass = camera.Cspad
    
        
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
        new = self.klass.from_basisgrid(bg)

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
        assert xyz2.shape == (32, 185, 388, 3)
        np.testing.assert_allclose(self.geom.xyz[0,0,0], xyz2[0])

        os.remove('ref_files/tmp_hdf.h5')