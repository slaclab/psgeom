
import numpy as np
import os
import unittest
import h5py

from psgeom import bin


def test_basic():
    
    pvs  = np.arange(10)
    data = np.random.randn(10)
    
    avg = bin.Averager(pvs, n_bins = 5)
    np.testing.assert_allclose(avg.bin_centers, np.arange(5) * 2)
    
    res = avg(data)
    np.testing.assert_allclose(res, (data[::2] + data[1::2]) / 2.0)

    return
    
    
def test_masking():

    pvs  = np.arange(10)
    data = np.random.randn(10)
    
    mask = np.ones(10, dtype=np.int)
    mask[::2] = 0
    
    avg = bin.Averager(pvs, mask, n_bins = 5)
    np.testing.assert_allclose(avg.bin_centers, np.arange(5) * 2)
    
    res = avg(data)
    np.testing.assert_allclose(res, data[1::2])
    
    return