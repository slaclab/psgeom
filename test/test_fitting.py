import numpy as np
import os
import unittest
import h5py

from psgeom import camera
from psgeom import fitting


class TestFitting(object):
    def test_basis_grid_interpolator(self):
        new_z = 0.75

        # load 3x geometires
        filenames = ['origin.geom', 'coffset05.geom', 'coffset10.geom']
        cameras = [camera.CompoundAreaCamera.from_crystfel_file('ref_files/distance_series/' + f) for f in filenames]
        motor_z = np.array([0.0, 0.5, 1.0])

        bgi = fitting.BasisGridInterpolator([g.to_basisgrid() for g in cameras], motor_z)
        prediction = bgi.predict( np.array([new_z]) )

        #print 'm:x (p/s/f)', bgi._coefficient_matrix[0,0::3]
        #print 'm:y (p/s/f)', bgi._coefficient_matrix[0,1::3]
        #print 'm:z (p/s/f)', bgi._coefficient_matrix[0,2::3]

        #print 'b:x (p/s/f)', bgi._coefficient_matrix[1,0::3]
        #print 'b:y (p/s/f)', bgi._coefficient_matrix[1,1::3]
        #print 'b:z (p/s/f)', bgi._coefficient_matrix[1,2::3]

        # the predicted bg should be the same as any other in x/y
        bg0 = cameras[0].to_basisgrid()
        xy_diff = np.sum(np.abs( prediction.xyz[...,:2] - bg0.xyz[...,:2] ))
        assert xy_diff < 1.0

        # and have a specific z-value
        #print prediction.xyz[...,2]
        bg2 = cameras[2].to_basisgrid()
        z_diff = (bg2.xyz[...,2] - bg0.xyz[...,2])
        z_expt = z_diff * 0.75 + bg0.xyz[...,2]
        #print z_expt
        z_diff = np.sum(np.abs(z_expt - prediction.xyz[...,2]))
        assert z_diff < 1.0

        return
