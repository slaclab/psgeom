# ---- moveable.py -------------------------------------------------------------

import numpy as np
import os
import unittest
import h5py

from psgeom import moveable


def test_translation_matrix_from_vector():
    
    x = np.random.randint(0,5,size=(3))    
    y = np.random.randint(0,5,size=(3))
    
    yp = np.ones(4)
    yp[:3] = y
    
    T = moveable._translation_matrix_from_vector(x)
    
    assert np.all(np.dot(T, yp)[:3] == x + y)
    assert np.dot(T, yp)[3] == 1.0
    
    
def test_rotation_matrix_from_angles():
    
    x = np.array([1.0, 0.0, 0.0]) # vector pointing at x
    
    Rz = moveable._rotation_matrix_from_angles(90.0, 0.0, 0.0)
    y = np.dot(Rz, x)
    np.testing.assert_array_almost_equal(y, np.array([0.0, 1.0, 0.0]),
                                         err_msg='Rz err')
    
    Rx = moveable._rotation_matrix_from_angles(0.0, 0.0, 90.0)
    z = np.dot(Rx, y)
    np.testing.assert_array_almost_equal(z, np.array([0.0, 0.0, 1.0]),
                                         err_msg='Rx err')
    
    Ry = moveable._rotation_matrix_from_angles(0.0, -90.0, 0.0)
    x = np.dot(Ry, z)
    np.testing.assert_array_almost_equal(x, x, err_msg='Ry err')
    

def test_angle_retrieval():

    for i in range(100):

        alpha = np.random.rand() * 180.0
        beta  = np.random.rand() * 180.0
        gamma = np.random.rand() * 180.0

        R = moveable._rotation_matrix_from_angles(gamma, beta, alpha)

        I = np.eye(3)
        Ip = np.dot(R, I)

        gp, bp, ap = moveable._angles_from_rotated_frame(Ip[:,0], Ip[:,1], Ip[:,2])

        #print np.array([gamma, beta, alpha]), np.array([gp, bp, ap])

        # test to make sure they rotate to the same thing
        R2 = moveable._rotation_matrix_from_angles(gp, bp, ap)

        assert np.sum( np.abs( R - R2 ) ) < 1e-12