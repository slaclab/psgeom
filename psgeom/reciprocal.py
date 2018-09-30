
import cPickle
import numpy as np

from psgeom import basisgrid
from psgeom.camera import arctan3

# ------------------------------------------------------------------------------
# FUNDAMENTAL CONSTANTS

h = 4.135677516e-15   # Planks constant | eV s
c = 299792458         # speed of light  | m / s

# ------------------------------------------------------------------------------


class PhotonEnergy(object):
    """
    Class that converts energies, wavelengths, frequencies, and wavenumbers.
    Each instance of this class can represent a light source.

    Attributes
    ----------
    self.energy      (eV)
    self.wavelength  (Angstroms)
    self.frequency   (Hz)
    self.wavenumber  (angular, inv. Angstroms)
    """

    def __init__(self, **kwargs):
        """
        An energy unit conversion class

        Parameters
        ----------
        **kwargs : dict
            Exactly one of the following, in the indicated units
            -- energy:     eV
            -- wavelength: angstroms
            -- frequency:  Hz
            -- wavenumber: inverse angstroms
        """

        # make sure we have only one argument
        if len(kwargs) != 1:
            raise KeyError('Expected exactly one argument, got %d' % (len(args)+1) )

        # no matter what gets provided, go straight to energy and then
        # convert to the rest from there
        for key in kwargs:

            if key == 'energy':
                self.energy = float(kwargs[key])

            elif key == 'wavenumber':
                self.wavenumber = float(kwargs[key])
                self.energy = self.wavenumber * h * c * 10.**10. / (2.0 * np.pi)

            elif key == 'wavelength':
                self.wavelength = float(kwargs[key])
                self.energy = h * c * 10.**10. / self.wavelength

            elif key == 'frequency':
                self.frequency = float(kwargs[key])
                self.energy = self.frequency * h

            else:
                raise ValueError('%s not a recognized kwarg' % key)

        # perform the rest of the conversions
        self.wavelength = h * c * 10.**10. / self.energy
        self.wavenumber = 2.0 * np.pi / self.wavelength
        self.frequency = self.energy / h

        # some aliases
        self.k = self.wavenumber


class Geometry(object):
    """
    Class that provides a plethora of geometric specifications for a detector
    setup.
    """

    def __init__(self, xyz, eV, beam_vector=None):
        """
        Geometry objects provide a handle for the many representations of
        detector geometry in scattering experiments, namely:

        -- real space
        -- real space in polar coordinates
        -- reciprocal space (q-space)
        -- reciprocal space in polar coordinates (q, theta, phi)

        Note the the origin is assumed to be the interaction site.

        Parameters
        ----------
        xyz : ndarray OR BasisGrid
            An a specification the (x,y,z) positions of each pixel. This can
            either be n x 3 array with the explicit positions of each pixel,
            or a BasisGrid object with a vectorized representation of the
            pixels. The latter yeilds higher performance, and is recommended.

        eV : float or reciprocal.PhotonEnergy
            The energy of the indicent beam, or a PhotonEnergy instance

        Optional Parameters
        -------------------
        beam_vector : float
            The 3-vector describing the beam direction. If `None`, then the
            beam is assumed to be purely in the z-direction.
        """

        if type(xyz) == np.ndarray:
            #print('xyz type: np.ndarray, initializing an explicit detector')

            self._pixels = xyz
            self._basis_grid = None
            self.num_pixels = xyz.shape[0]
            self._xyz_type = 'explicit'

        elif type(xyz) == basisgrid.BasisGrid:
            #print('xyz type: BasisGrid, initializing an implicit detector')

            self._pixels = None
            self._basis_grid = xyz
            self.num_pixels = self._basis_grid.num_pixels
            self._xyz_type = 'implicit'

        else:
            raise TypeError("`xyz` type must be one of {'np.ndarray', "
                            "'thor.xray.BasisGrid'}")


        # parse wavenumber
        if isinstance(eV, PhotonEnergy):
            self.energy = eV
        elif type(eV) in [float, np.float64, np.float32]:
            self.energy = PhotonEnergy(energy=eV)
        else:
            raise TypeError('`eV` must be a float or '
                            'psgeom.reciprocal.PhotonEnergy')
        self.k = self.energy.wavenumber


        # parse beam_vector -- is guarenteed to be a unit vector
        if beam_vector is not None:
            if beam_vector.shape == (3,):
                self.beam_vector = self._unit_vector(beam_vector)
            else:
                raise ValueError('`beam_vector` must be a 3-vector')
        else:
            self.beam_vector = np.array([0.0, 0.0, 1.0])

        return


    def implicit_to_explicit(self):
        """
        Convert an implicit Geometry to an explicit one (where the xyz pixels
        are stored in memory).
        """
        if not self.xyz_type == 'implicit':
            raise Exception('Geometry must have xyz_type implicit for conversion.')
        self._pixels = self.xyz
        self._xyz_type = 'explicit'
        return


    @property
    def xyz_type(self):
        return self._xyz_type


    @property
    def xyz(self):
        if self.xyz_type == 'explicit':
            return self._pixels
        elif self.xyz_type == 'implicit':
            return self._basis_grid.to_explicit()


    @property
    def real(self):
        return self.xyz.copy()


    @property
    def polar(self):
        return self._real_to_polar(self.real)


    @property
    def reciprocal(self):
        return self._real_to_reciprocal(self.real)


    @property
    def recpolar(self):
        a = self._real_to_recpolar(self.real)
        # convention: theta is angle of q-vec with plane normal to beam
        a[:,1] = self.polar[:,1] / 2.0
        return a


    @property
    def q_max(self):
        """
        Returns the maximum value of |q| the detector measures
        """

        if self.xyz_type == 'explicit':
            q_max = np.max(self.recpolar[:,0])

        elif self.xyz_type == 'implicit':
            q_max = 0.0
            for i in range(self._basis_grid.num_grids):
                c  = self._basis_grid.get_grid_corners(i)
                qc = self._real_to_recpolar(c)
                q_max = max([q_max, float(np.max(qc[:,0]))])

        return q_max


    def evaluate_qmag(self, xyz):
        """
        Given the positions of pixels `xyz`, compute the corresponding |q|
        value for each.

        Parameters
        ----------
        qxyz : ndarray, float
            The array of pixels (shape : N x 3)

        Returns
        -------
        qmag : ndarray, float
            The array of q-magnitudes, len N.
        """
        thetas = self._evaluate_theta(xyz)
        qmag = 2.0 * self.k * np.sin(thetas/2.0)
        return qmag


    def _evaluate_theta(self, xyz):
        """
        Given the positions of pixels `xyz`, compute the corresponding
        scattering angle theta for each.

        Parameters
        ----------
        xyz : ndarray, float
            The array of pixels (shape : N x 3)

        Returns
        -------
        thetas : ndarray, float
            The scattering angles for each pixel
        """
        u_xyz  = self._unit_vector(xyz)
        thetas = np.arccos(np.dot( u_xyz, self.beam_vector ))
        return thetas


    def _real_to_polar(self, xyz):
        """
        Convert the real-space representation to polar coordinates.
        """
        polar = self._to_polar(xyz)
        return polar


    def _real_to_reciprocal(self, xyz):
        """
        Convert the real-space to reciprocal-space in cartesian form.
        """

        assert len(xyz.shape) == 2
        assert xyz.shape[1] == 3

        # generate unit vectors in the pixel direction, origin at sample
        S = self._unit_vector(xyz)
        q = self.k * (S - self.beam_vector)

        return q


    def _real_to_recpolar(self, xyz):
        """
        Convert the real-space to reciprocal-space in polar form, that is
        (|q|, theta , phi).
        """
        reciprocal_polar = self._to_polar( self._real_to_reciprocal(xyz) )
        return reciprocal_polar


    @staticmethod
    def _norm(vector):
        """
        Compute the norm of an n x m array of vectors, where m is the dimension.
        """
        if len(vector.shape) == 2:
            assert vector.shape[1] == 3
            norm = np.sqrt( np.sum( np.power(vector, 2), axis=1 ) )
        elif len(vector.shape) == 1:
            assert vector.shape[0] == 3
            norm = np.sqrt( np.sum( np.power(vector, 2) ) )
        else:
            raise ValueError('Shape of vector wrong')
        return norm


    def _unit_vector(self, vector):
        """
        Returns a unit-norm version of `vector`.

        Parameters
        ----------
        vector : ndarray, float
            An n x m vector of floats, where m is assumed to be the dimension
            of the space.

        Returns
        -------
        unit_vectors : ndarray,float
            An n x m vector, same as before, but now of unit length
        """

        norm = self._norm(vector)

        if len(vector.shape) == 1:
            unit_vectors = vector / norm

        elif len(vector.shape) == 2:
            unit_vectors = np.zeros( vector.shape )
            for i in range(vector.shape[0]):
                unit_vectors[i,:] = vector[i,:] / norm[i]

        else:
            raise ValueError('invalid shape for `vector`: %s' % str(vector.shape))

        return unit_vectors


    def _to_polar(self, vector):
        """
        Converts n m-dimensional `vector`s to polar coordinates. By polar
        coordinates, I mean the cannonical physicist's (r, theta, phi), no
        2-theta business.

        We take, as convention, the 'z' direction to be along self.beam_vector
        """

        polar = np.zeros( vector.shape )

        # note the below is a little modified from the standard, to take into
        # account the fact that the beam may not be only in the z direction

        polar[:,0] = self._norm(vector)
        polar[:,1] = np.arccos( np.dot(vector, self.beam_vector) / \
                                (polar[:,0]+1e-16) )           # cos^{-1}(z.x/r)
        polar[:,2] = arctan3(vector[:,1] - self.beam_vector[1],
                             vector[:,0] - self.beam_vector[0])   # y first!

        return polar


    def compute_intersections(self, q_vectors, grid_index, run_checks=True):
        """
        Compute the points i=(x,y,z) where the scattering vectors described by
        `q_vectors` intersect the detector.

        Parameters
        ----------
        q_vectors : np.ndarray
            An N x 3 array representing q-vectors in cartesian q-space.

        grid_index : int
            The index of the grid array to intersect

        Optional Parameters
        -------------------
        run_checks: bool
            Whether to run some good sanity checks, at small computational cost.

        Returns
        -------
        pix_n : ndarray, float
            The coefficients of the position of each intersection in terms of
            the basis grids s/f vectors.

        intersect: ndarray, bool
            A boolean array of which of `q_vectors` intersect with the grid
            plane. First column is slow scan vector (s),  second is fast (f).

        References
        ----------
        .[1] http://en.wikipedia.org/wiki/Line-plane_intersection
        """

        if not self.xyz_type == 'implicit':
            raise RuntimeError('intersections can only be computed for implicit'
                               ' detectors')

        # compute the scattering vectors corresponding to q_vectors
        S = (q_vectors / self.k) + self.beam_vector

        # compute intersections
        p, s, f, shape = self._basis_grid.get_grid(grid_index)
        n = self._unit_vector( np.cross(s, f) )
        i = (np.dot(p, n) / np.dot(S, n))[:,None] * S

        # convert to pixel units by solving for the coefficients of proj
        A = np.array([s,f]).T
        pix_n, resid, rank, sigma = np.linalg.lstsq( A, (i-p).T )

        if run_checks:
            err = np.sum( np.abs((i-p) - np.transpose( np.dot(A, pix_n) )) )
            if err > 1e-6:
                raise RuntimeError('Error in computing where scattering vectors '
                                   'intersect with detector. Intersect not reproduced'
                                   ' (err: %f per pixel)' % (err / i.shape[0],) )

            if not np.sum(resid) < 1e-6:
                raise RuntimeError('Error in basis grid (residuals of point '
                                   'placement too large). Perhaps fast/slow '
                                   'vectors describing basis grid are linearly '
                                   'dependant?')

        pix_n = pix_n.T
        assert pix_n.shape[1] == 2 # s/f

        # see if the intersection in the plane is on the detector grid
        intersect = (pix_n[:,0] >= 0.0) * (pix_n[:,0] <= float(shape[0]-1)) *\
                    (pix_n[:,1] >= 0.0) * (pix_n[:,1] <= float(shape[1]-1))

        #print('%.3f %% of pixels intersected by grid %d' % \
        #    ( (float(np.sum(intersect)) / float(np.product(intersect.shape)) * 100.0),
        #    grid_index) )

        return pix_n[intersect], intersect


    @classmethod
    def generic(cls, spacing=1.00, lim=100.0, eV=10000.0, l=50.0,
                force_explicit=False):
        """
        Generates a simple grid detector that can be used for testing
        (factory function).

        Optional Parameters
        -------------------
        spacing : float
            The real-space grid spacing
        lim : float
            The upper and lower limits of the grid
        eV : float
            Energy of the beam (in eV)
        l : float
            The path length from the sample to the detector, in the same units
            as the detector dimensions.
        force_explicit : bool
            Forces the detector to be xyz_type explicit. Mostly for debugging.
            Recommend keeping `False`.

        Returns
        -------
        geometry : reciprocal.Geometry
            An instance of the detector that meets the specifications of the
            parameters
        """

        if not force_explicit:

            p = np.array([-lim, -lim, l])   # corner position
            f = np.array([0.0, spacing, 0.0]) # slow scan is x
            s = np.array([spacing, 0.0, 0.0]) # fast scan is y

            dim = int(2*(lim / spacing) + 1)
            shape = (dim, dim)

            basis = basisgrid.BasisGrid()
            basis.add_grid(p, s, f, shape)

            geo = cls(basis, eV)

        else:
            x = np.arange(-lim, lim+spacing, spacing)
            xx, yy = np.meshgrid(x, x)

            xyz = np.zeros((len(x)**2, 3))
            xyz[:,0] = yy.flatten() # fast scan is y
            xyz[:,1] = xx.flatten() # slow scan is x
            xyz[:,2] = l

            geo = cls(xyz, eV)

        return geo


    def _to_serial(self):
        """ serialize the object to an array """
        s = np.array( cPickle.dumps(self) )
        s.shape=(1,) # a bit nasty...
        return s


    @classmethod
    def _from_serial(self, serialized):
        """ recover a Geometry object from a serialized array """
        if serialized.shape == (1,):
            serialized = serialized[0]
        d = cPickle.loads( str(serialized) )
        return d


