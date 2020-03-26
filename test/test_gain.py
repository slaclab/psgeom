# --- gain.py -----------------------------------------------------------------

import numpy as np
import os
import unittest
import h5py

from psgeom import gain

class TestGain:
    
    def test_cheetah_roundtrip(self):
        
        r = np.random.randint(2, size=(32,185,388))
        gain.write_cheetah('../ref_filestmpgain.h5', r)
        r2 = gain.load_cheetah('../ref_filestmpgain.h5')
        
        e = np.sum(np.abs( r - r2 ))
        assert e == 0, e
        os.remove('../ref_filestmpgain.h5')
        
        return
        
    def test_DAQ_roundtrip(self):

        # DAQ map should be comprised of 1's and 7.2's
        r = np.random.randint(2, size=(32,185,388)) * 6.2 + 1
        gain.write_daq('../ref_filestmpgain.txt', r)
        r2 = gain.load_daq('../ref_filestmpgain.txt')
        
        e = np.sum(np.abs( r - r2 ))
        np.testing.assert_allclose(r, r2)
        os.remove('../ref_filestmpgain.txt')

        return
    
    def test_translate_consistency(self):
        
        cht = gain.load_cheetah('ref_files/200px-gainmap.h5')
        daq = gain.load_daq('ref_files/200px-gainmap.txt', gain_ratio=7.2)
        
        assert np.sum(np.abs(cht - daq)) == 0.0
        
        return