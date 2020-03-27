# --- metrology.py ------------------------------------------------------------

import numpy as np
import os
import unittest
import h5py

from psgeom import metrology
from psgeom import camera

def test_load_metrology():

    # To do: these are just smoke tests, can we make them
    #        real tests?

    x  = metrology.load_to_basisgrid('ref_files/cspad/refgeom_metrology.txt')
    assert x.xyz.shape == (2296960, 3)

    x2 = camera.Cspad.from_metrology_file('ref_files/cspad/refgeom_metrology.txt')
    assert x2.xyz.shape == (4, 8, 185, 388, 3)
    #np.testing.assert_allclose(x.xyz, x2.xyz)

    return