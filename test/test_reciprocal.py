# --- reciprocal.py -----------------------------------------------------------
#   + basisgrid.py
#   + fitting.py

import numpy as np
import os
import unittest
import h5py

from psgeom import basisgrid
from psgeom import reciprocal
from psgeom import camera
from psgeom import fitting


def test_bg_as_array():
    # prob not necessary
    geom = camera.Cspad.from_psana_file('ref_files/cspad/refgeom_psana.data')
    bg = geom.to_basisgrid()
    assert bg.as_array().shape == (64, 11)

def test_bg_from_array():
    geom = camera.Cspad.from_psana_file('ref_files/cspad/refgeom_psana.data')
    bg = geom.to_basisgrid()
    bg2 = basisgrid.BasisGrid.from_array( bg.as_array() )
    assert np.all( bg.to_explicit() == bg2.to_explicit() )
    


class TestPhotonEnergy(object):

    def test_unit_convs(self):
        b = reciprocal.PhotonEnergy(energy=1000.0) # 1 keV
        np.testing.assert_allclose(b.wavelength, 12.398, rtol=1e-3)
        np.testing.assert_allclose(b.frequency, 2.4190e17, rtol=1e-3)
        np.testing.assert_allclose(b.wavenumber, (2.0 * np.pi)/12.398, rtol=1e-3)


class TestGeometry(object):

    def setup(self):
        self.spacing   = 0.05
        self.lim       = 10.0
        self.energy    = reciprocal.PhotonEnergy(wavenumber=0.7293)
        self.l         = 50.0
        self.d = reciprocal.Geometry.generic(spacing = self.spacing,
                                             lim = self.lim,
                                             eV = self.energy,
                                             l = self.l)

    def test_implicit_to_explicit(self):
        xyz_imp = self.d.real
        self.d.implicit_to_explicit()
        np.testing.assert_array_almost_equal(xyz_imp, self.d.real)

    def test_evaluate_qmag(self):
        # doubles as a test for _evaluate_theta
        x = np.zeros((5, 3))
        x[:,0] = np.random.randn(5)
        x[:,2] = self.l

        S = x.copy()
        S = S / np.sqrt( np.sum( np.power(S, 2), axis=1 ) )[:,None]
        S -= self.d.beam_vector

        b = reciprocal.PhotonEnergy(wavenumber=0.7293)
        qref = b.k * np.sqrt( np.sum( np.power(S, 2), axis=1 ) )

        qmag = self.d.evaluate_qmag(x)
        np.testing.assert_allclose(qref, qmag)

    def test_recpolar_n_reciprocal(self):
        q1 = np.sqrt( np.sum( np.power(self.d.reciprocal,2), axis=1) )
        q2 = self.d.recpolar[:,0]
        np.testing.assert_array_almost_equal(q1, q2)

    def test_polar_space(self):

        # this is the "generic" geometry in real space
        x = np.arange(-self.lim, self.lim+self.spacing, self.spacing)
        xx, yy = np.meshgrid(x, x)

        # one slice along the horizontal direction in real space
        r     = self.d.polar[:,0]
        theta = self.d.polar[:,1]
        phi   = self.d.polar[:,2]

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        np.testing.assert_array_almost_equal(yy.flatten(), x)
        np.testing.assert_array_almost_equal(xx.flatten(), y)

    def test_reciprocal_space(self):
        qx = self.d.reciprocal[:,0]
        qy = self.d.reciprocal[:,1]

        Shat    = self.d._unit_vector(self.d.real)
        Sx_unit = Shat[:,0]
        Sy_unit = Shat[:,1]

        np.testing.assert_array_almost_equal(qx/self.d.k, Sx_unit)
        np.testing.assert_array_almost_equal(qy/self.d.k, Sy_unit)

    def test_recpolar_space(self):

        # build a reference conversion, using a different geometrical calc
        ref1 = np.zeros(self.d.xyz.shape)
        hd = np.sqrt( np.power(self.d.xyz[:,0], 2) + np.power(self.d.xyz[:,1], 2) )

        # |q| = k*sqrt{ 2 - 2 cos(theta) }
        ref1[:,0] = self.d.k * np.sqrt( 2.0 - 2.0 * np.cos(self.d.polar[:,1]) )

        # q_theta = theta / 2 (one-theta convention)
        ref1[:,1] = self.d.polar[:,1] / 2.0 # not working atm

        # q_phi is the same as polar
        ref1[:,2] = self.d.polar[:,2].copy()

        np.testing.assert_array_almost_equal(ref1[:,0], self.d.recpolar[:,0], err_msg='|q|')
        np.testing.assert_array_almost_equal(ref1[:,1], self.d.recpolar[:,1], err_msg='theta')
        np.testing.assert_array_almost_equal(ref1[:,2], self.d.recpolar[:,2], err_msg='phi')

    def test_compute_intersect(self):

        # build a simple grid and turn it into a geometry
        bg = basisgrid.BasisGrid()
        p = np.array([0.0, 0.0, 1.0])
        s = np.array([1.0, 0.0, 0.0])
        f = np.array([0.0, 1.0, 0.0])
        shape = (10, 10)
        bg.add_grid(p, s, f, shape)
        d = reciprocal.Geometry(bg, 2.0*np.pi/1.4)

        # compute a set of q-vectors corresponding to a slightly offset grid
        xyz_grid = bg.to_explicit()
        xyz_off = xyz_grid.copy()
        xyz_off[:,0] += 0.5
        xyz_off[:,1] += 0.5
        q_vectors = d._real_to_reciprocal(xyz_off)

        # b/c s/f vectors are unit vectors, where they intersect s/f is simply
        # their coordinates. The last row and column will miss, however
        intersect_ref = np.logical_and( (xyz_off[:,0] <= 9.0),
                                        (xyz_off[:,1] <= 9.0) )

        pix_ref = xyz_off[intersect_ref,:2]

        # compute the intersection from code
        pix, intersect = d.compute_intersections(q_vectors, 0) # 0 --> grid_index
        print(pix, intersect)

        np.testing.assert_array_almost_equal(intersect_ref, intersect)
        np.testing.assert_array_almost_equal(pix_ref, pix)

    def test_serialization(self):
        s = self.d._to_serial()
        d2 = reciprocal.Geometry._from_serial(s)
        np.testing.assert_array_almost_equal(d2.xyz, self.d.xyz)

    def test_q_max(self):
        ref_q_max = np.max(self.d.recpolar[:,0])
        np.testing.assert_almost_equal(self.d.q_max, ref_q_max, decimal=2)
        
        
class TestFitting(object):
    def test_basis_grid_interpolator(self):
        new_z = 0.75

        # load 3x geometires
        filenames = ['origin.geom', 'coffset05.geom', 'coffset10.geom']
        cameras = [camera.Cspad.from_crystfel_file('ref_files/distance_series/' + f) for f in filenames]
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