# ---- translate.py ------------------------------------------------------------
    
import numpy as np
import os
import unittest
import h5py

from psgeom import camera
from psgeom import translate 
from psgeom import sensors   
    
PIXEL_TOLERANCE_um = 10.0


class TestTranslate(object):
    """
    Will test read/write from external software packages
    """
    
    def setup(self):
        self.cd = camera.CompoundCamera.from_psana_file('ref_files/cspad/refgeom_psana.data')
        self.cspad = camera.Cspad.from_psana_file('ref_files/cspad/refgeom_psana.data')
                
        
    def test_psf_text(self):
        
        raise unittest.SkipTest
        
        self.cd.to_text_file('ref_files/cd_psf.txt')
        self.cspad.to_text_file('ref_files/cspad_psf.txt')
        
        # todo : load & ensure consistent
        cd2 = camera.CompoundCamera.from_text_file('ref_files/cd_psf.txt')
        cspad2 = camera.CompoundCamera.from_text_file('ref_files/cspad_psf.txt')
        
        np.testing.assert_allclose(self.cd.xyz, cd2.xyz)
        np.testing.assert_allclose(self.cspad.xyz, cspad2.xyz)
        
        os.remove('ref_files/cd_psf.txt')
        os.remove('ref_files/cspad_psf.txt')

        
    def test_cheetah_roundtrip(self):
        
        self.cspad.to_cheetah_file('ref_files/tmp_cheetah_geom.h5')
        cspad2 = camera.Cspad.from_cheetah_file('ref_files/tmp_cheetah_geom.h5')
                
        np.testing.assert_allclose( np.squeeze(self.cspad.xyz),
                                    np.squeeze(cspad2.xyz),
                                    err_msg='round trip fail',
                                    atol=PIXEL_TOLERANCE_um)
        os.remove('ref_files/tmp_cheetah_geom.h5')
        
        
    def test_cheetah_values(self):
        
        # the cheetah pixels are just the xyz, so we can load them and use them
        # to compare to our implementation
        
        f = h5py.File('ref_files/cspad/refgeom_cheetah.h5', 'r')
        
        x = -1.0 * translate._cheetah_to_twobyones( np.array(f['x']) * 1000000.0 )
        y = translate._cheetah_to_twobyones( np.array(f['y']) * 1000000.0 )
        z = translate._cheetah_to_twobyones( np.array(f['z']) )
        
        f.close()
        
        ref_xyz = np.rollaxis( np.array([x,y,z]), 0, 4) # to shape (32,...,3)
        tst = camera.Cspad.from_cheetah_file('ref_files/cspad/refgeom_cheetah.h5')
        tst_xyz = tst.xyz.reshape(32,185,388,3)
                
        np.testing.assert_allclose(tst_xyz[:,0,0,:],
                                   ref_xyz[:,0,0,:],
                                   atol=PIXEL_TOLERANCE_um,
                                   err_msg='panel 1st pixels off')
        
        
                                   
        # cheetah does not deal correctly with the large center pixels, so
        # we test around that
        assert np.max( np.abs(tst_xyz - ref_xyz) ) < 500.0
        np.testing.assert_allclose(tst_xyz[:,:,:,:], ref_xyz[:,:,:,:],
                                   atol=PIXEL_TOLERANCE_um,
                                   err_msg='panels off in general')
        
        
    def test_cheetah_crystfel_consistency(self):
        
        # the following two files were generated using the "make_pixelmap"
        # program by TAW, the cheetah file being generated from the CrystFEL
        # conversion done by TJL 7/10/15
        
        # very close but ~1/2 pixel off. Does cheetah plot pixel corners or
        # centers? -- TJL 7/8/15
        
        cheetah  = camera.Cspad.from_cheetah_file('ref_files/cspad/refgeom_cheetah.h5')
        crystfel = camera.Cspad.from_crystfel_file('ref_files/cspad/refgeom_crystfel.geom')
        
        print(np.squeeze(cheetah.xyz) - np.squeeze(crystfel.xyz))
        
        # cheetah z-values are screwed up, so ignore those
        np.testing.assert_allclose(np.squeeze(cheetah.xyz)[...,:2],
                                   np.squeeze(crystfel.xyz)[...,:2],
                                   atol=PIXEL_TOLERANCE_um)
       
        
    def test_crystfel_roundtrip(self):
        
        # first test with no weird z stuff
        
        self.cspad.to_crystfel_file('ref_files/tmp_crystfel.geom')
        cd2 = camera.Cspad.from_crystfel_file('ref_files/tmp_crystfel.geom')
        
        print('HERE1', self.cspad.xyz.shape, cd2.xyz.shape)
        
        # be sure error is less than 1 micron in x/y, 0.2 mm in z
        assert np.max(np.abs( np.squeeze(self.cspad.xyz[...,:2]) - np.squeeze(cd2.xyz[...,:2]) )) < 1.0
        assert np.max(np.abs( np.squeeze(self.cspad.xyz) - np.squeeze(cd2.xyz) )) < 200.0

        print('HERE2')        
        
        # CrystFEL's geometry assumes all panels are orthogonal to the beam. To
        # do a fair comparison, therefore, we set all the z-rotations to zero
        # note: we remove rotations around the x- & y-axes to get rid of z
        
        for c in self.cspad.children:
            c._rotation_angles[1:] = 0.0
            for q in c.children:
                q._rotation_angles[1:] = 0.0
                for t in q.children:
                    t._rotation_angles[1:] = 0.0
        self.cspad._rotation_angles[1:] = 0.0
        
        for l in self.cspad.leaves:
            assert [v[2] == 0.0 for v in l.psf[1:]]
        
        
        self.cspad.to_crystfel_file('ref_files/tmp_crystfel.geom')
        cd2 = camera.Cspad.from_crystfel_file('ref_files/tmp_crystfel.geom')
        
        assert np.max(np.abs( np.squeeze(self.cspad.xyz) - np.squeeze(cd2.xyz) )) < 1.0
        
        
        # finally just check all pixels are straight up reasonable
        np.testing.assert_allclose(np.squeeze(self.cd.xyz),
                                   np.squeeze(cd2.xyz),
                                   err_msg='round trip fail',
                                   atol=PIXEL_TOLERANCE_um,
                                   rtol=1e-3)
                                           
        os.remove('ref_files/tmp_crystfel.geom')
        
        
    def test_crystfel_coffset(self):
        
        geom = camera.Cspad.from_psana_file('ref_files/cspad/refgeom_psana.data')
        geom.to_crystfel_file('ref_files/temp.geom', coffset=2.0)

        geom2 = camera.Cspad.from_crystfel_file('ref_files/temp.geom')

        assert np.all(np.abs(geom2.xyz[...,2] - 2e6) < 0.01)
        np.testing.assert_allclose(np.squeeze(geom.xyz[...,:2]),
                                   np.squeeze(geom2.xyz[...,:2]),
                                   err_msg='round trip fail',
                                   atol=PIXEL_TOLERANCE_um,
                                   rtol=1e-3)



        os.remove('ref_files/temp.geom')


    def test_crystefel_nonstandard(self):
        # (1) panel names in crystfel may not begin with "p"
        # (2) ss and fs fields may have only x or y, and assume the other is 0.0
        
        geom = camera.CompoundAreaCamera.from_crystfel_file('ref_files/pnccd/pnccd2.geom')
        print(geom.xyz.shape)
        assert geom.xyz.shape == (2, 1024, 512, 3), geom.xyz.shape
        geom.to_crystfel_file('ref_files/temp.geom', coffset=2.0)
        os.remove('ref_files/temp.geom')

    
    def test_2x2_cheetah_roundtrip(self):
        
        cheetah  = camera.Cspad.from_cheetah_file('ref_files/cspad-2x2/cspad-2x2-approx1.h5')
        cheetah.to_psana_file('ref_files/tmp_psana.data')
        cheetah2 = camera.Cspad.from_psana_file('ref_files/tmp_psana.data')
        
        np.testing.assert_allclose(np.squeeze(cheetah.xyz),
                                   np.squeeze(cheetah2.xyz),
                                   atol=PIXEL_TOLERANCE_um)
                                   
        os.remove('ref_files/tmp_psana.data')
        
        
    def test_2x2_crystfel_roundtrip(self):
        
        crystfel  = camera.Cspad.from_crystfel_file('ref_files/cspad-2x2/cspad-2x2-approx1.geom')
        crystfel.to_psana_file('ref_files/tmp_psana.data')
        crystfel2 = camera.Cspad.from_psana_file('ref_files/tmp_psana.data')
        
        np.testing.assert_allclose(np.squeeze(crystfel.xyz),
                                   np.squeeze(crystfel2.xyz),
                                   atol=PIXEL_TOLERANCE_um)
                                   
        os.remove('ref_files/tmp_psana.data')
        
        
    def test_2x2_consistency(self):
        
        # TJL: I have not verified that these are in fact the same geom
        #      they were provided by users
        crystfel = camera.Cspad.from_crystfel_file('ref_files/cspad-2x2/cspad-2x2-approx1.geom')
        cheetah  = camera.Cspad.from_cheetah_file('ref_files/cspad-2x2/cspad-2x2-approx1.h5')
        
        # compare only x/y, not z
        np.testing.assert_allclose(np.squeeze(crystfel.xyz)[...,:2],
                                   np.squeeze(cheetah.xyz)[...,:2],
                                   atol=100.0)


    def test_dials_load(self):
        # TODO make this a real test
        obj = camera.Cspad()
        a = translate.load_dials(obj, 'ref_files/cspad/refgeom_dials.json')
        assert a.xyz.shape == (4, 8, 185, 388, 3)

        obj2 = camera.CompoundAreaCamera()
        x = translate.load_dials(obj2, 'ref_files/cspad/refgeom_dials2.json')
        assert x.xyz.shape == (8, 512, 1024, 3)


    def test_jungfrau1M_ref_file(self):
        
        obj = camera.CompoundAreaCamera.from_crystfel_file('ref_files/jungfrau/jungfrau-1M-pan.geom',
                                                           element_type=sensors.JungfrauSegment)
        assert obj.xyz.shape == (2, 512, 1024, 3)
        
        bg = obj.to_basisgrid()
        assert bg.num_grids == 16
        
        # roundtrip...
        obj.to_crystfel_file('tmp.geom')
        obj2 = camera.CompoundAreaCamera.from_crystfel_file('tmp.geom', element_type=sensors.JungfrauSegment)
        np.testing.assert_allclose(obj.xyz, obj2.xyz, atol=0.01)
        
        os.remove('tmp.geom')
        
        