# --- scripts ------------------------------------------------------------------

import os
import h5py
import numpy as np

from psgeom import camera

class TestGeoconv(object):
    
    def setup(self):
        self.cd = camera.CompoundCamera.from_psana_file('ref_files/refgeom_psana.data')
        self.cspad = camera.Cspad.from_psana_file('ref_files/refgeom_psana.data')
    
    def test_psana(self):
        os.system("geoconv --cspad -f psana ./ref_files/refgeom_psana.data tmp.data")
        cspad2 = camera.Cspad.from_psana_file('tmp.data')
        assert np.max(np.abs( np.squeeze(self.cspad.xyz[...,:2]) - np.squeeze(cspad2.xyz[...,:2]) )) < 1.0
        os.remove("tmp.data")
                
    def test_crystfel(self):
        os.system("geoconv --cspad -f crystfel ./ref_files/refgeom_psana.data tmp.geom")
        cspad2 = camera.Cspad.from_crystfel_file('tmp.geom')
        assert np.max(np.abs( np.squeeze(self.cspad.xyz[...,:2]) - np.squeeze(cspad2.xyz[...,:2]) )) < 1.0
        os.remove("tmp.geom")

    def test_cheetah(self):
        os.system("geoconv --cspad -f cheetah ./ref_files/refgeom_psana.data tmp.h5")
        cspad2 = camera.Cspad.from_cheetah_file('tmp.h5')
        assert np.max(np.abs( np.squeeze(self.cspad.xyz[...,:2]) - np.squeeze(cspad2.xyz[...,:2]) )) < 1.0
        os.remove("tmp.h5")
        
    def test_hdf5(self):
        os.system("geoconv --cspad -f hdf5 ./ref_files/refgeom_psana.data tmp.hdf5")
        with h5py.File('tmp.hdf5') as f:
            cspad2 = f['xyz']
            assert np.max(np.abs( np.squeeze(self.cd.xyz[...,:2]).reshape(32,185,388,2) \
               - np.squeeze(cspad2[...,:2]) )) < 1.0
        os.remove("tmp.hdf5")
        