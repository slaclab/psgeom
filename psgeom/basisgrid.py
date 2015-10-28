
"""
basisgrid.py
"""

import numpy as np

class BasisGrid(object):
    """
    A class representing a set of rectangular grids in space -- specifically,
    x-ray scattering detectors. Does not contain all the metadata associated
    with a full-fledged Detector class (e.g. the wavelength, etc).

    Note that the geometry below is definied in "slow" and "fast" scan
    dimensions. These are simply the two dimensions that define the plane
    a single rectangular pixel grid lives in. They may also be called the y and
    x dimensions without any loss of generality.

    The convention here -- and in all of ODIN -- is one of Row-Major ordering,
    which is consistent with C/python. This means that y is the slow dim, x is
    the fast dim, and when ordering these two they will appear as (slow, fast).

    Note on units: units are arbitrary -- all the units must be the same for
    this to work. We don't keep track of units here.

    The class is a set of rectangular grids, with each grid defined by four
    quantities:

        -- p vector : DEFINES A GRIDS POSITION IN SPACE.
                      The vector between a chosen origin (possibly interaction
                      site) and the corner of the grid that is smallest in both
                      slow and fast (x/y) dimensions of the coordinate system.
                      Usually this will be the "bottom left" corner, but due to
                      the generality of the coordinates used, this is not
                      necessarily true.

        -- s/f vect : DEFINES A GRIDS ORIENTATION IN SPACE
                      Vectors pointing along the slow/fast-scan direction,
                      respectively. These define the plane containing the pixels.
                      The magnitudes of these vectors defines the size of the
                      pixel in that dimension.

        -- shape    : DEFINES GRID DIMENSIONS
                      The number of pixels in the fast/slow direction. Ints.
    """


    def __init__(self, list_of_grids=[]):
        """
        Initialize a BasisGrid object.

        Parameters
        ----------
        list_of_grids : list
            A list of tuples of the form  (p, s, f, shape). See the doc
            for the `add_grid` method on this class for more information. May
            be an empty list (default) in which case a GridList with no pixels
            is created.

        See Also
        --------
        add_grid
        add_grid_using_center
        """

        if not type(list_of_grids) == list:
            raise TypeError('`list_of_grids` must be a list')

        self._num_grids = 0
        self._ps        = [] # p-vectors
        self._ss        = [] # slow-scan vectors
        self._fs        = [] # fast-scan vectors
        self._shapes    = [] # shapes

        if len(list_of_grids) > 0:
            for grid in list_of_grids:
                self.add_grid(*grid)

        return
        

    def _check_valid_basis(self, p, s, f, shape):
        """
        Check to make sure that all the inputs look good.
        """

        if not (p.shape == (3,)) and (s.shape == (3,)) and (f.shape == (3,)):
            raise ValueError('`p`, `s`, `f` must be 3-vectors')

        if not (len(shape) == 2):
            raise ValueError('`shape` must be len 2')

        return


    def _assert_list_sizes(self):
        """
        A simple sanity check
        """
        assert len(self._ps)     == self.num_grids
        assert len(self._ss)     == self.num_grids
        assert len(self._fs)     == self.num_grids
        assert len(self._shapes) == self.num_grids
        return


    @property
    def num_pixels(self):
        """
        Return the total number of pixels in the BasisGrid.
        """
        n = np.sum([np.product(self._shapes[i]) for i in range(self.num_grids)])
        return int(n)


    @property
    def num_grids(self):
        return self._num_grids


    def add_grid(self, p, s, f, shape):
        """
        Add a grid (detector array) to the basis representation.

        Parameters
        ----------
        p : np.ndarray, float
            3-vector from the origin to the pixel on the grid with
            smallest coordinate in all dimensions.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.

        See Also
        --------
        add_grid_using_center
        """
        self._check_valid_basis(p, s, f, shape)
        self._ps.append(p)
        self._ss.append(s)
        self._fs.append(f)
        self._shapes.append(shape)
        self._num_grids += 1
        self._assert_list_sizes()
        return


    def add_grid_using_center(self, p_center, s, f, shape):
        """
        Add a grid (detector array) to the basis representation. Here, though,
        the p-vector points to the center of the array instead of the slow/fast
        smallest corner.

        Parameters
        ----------
        p_center : np.ndarray, float
            3-vector from the origin to the center of the grid.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.
        """

        p_center = np.array(p_center)
        if not p_center.shape == (3,):
            raise ValueError('`p_center` must have shape (3,)')

        # just compute where `p` is then add the grid as usual
        x = (np.array(shape) - 1)
        center_correction =  ((x[0] * s) + (x[1] * f)) / 2.
        p  = p_center.copy()
        p -= center_correction

        self.add_grid(p, s, f, shape)

        return


    def get_grid(self, grid_number):
        """
        Return a grid for grid `grid_number`.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get.

        Returns
        -------
        p_center : np.ndarray, float
            3-vector from the origin to the center of the grid.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.
        """

        if grid_number >= self.num_grids:
            raise ValueError('Only %d grids in object, you asked for the %d-th'
                             ' (zero indexed)' % (self.num_grids, grid_number))

        grid_tuple = (self._ps[grid_number], self._ss[grid_number],
                      self._fs[grid_number], self._shapes[grid_number])

        return grid_tuple


    def get_grid_corners(self, grid_number):
        """
        Return the positions of the four corners of a grid.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get the corners of.

        Returns
        -------
        corners : np.ndarray, float
            A 4 x 3 array, where the first dim represents the four corners, and
            the second is x/y/z. Note one corner is always just the `p` vector.
        """

        if grid_number >= self.num_grids:
            raise ValueError('Only %d grids in object, you asked for the %d-th'
                             ' (zero indexed)' % (self.num_grids, grid_number))

        # compute the lengths of the parallelogram sides
        s_side = self._ss[grid_number] * float(self._shapes[grid_number][0])
        f_side = self._fs[grid_number] * float(self._shapes[grid_number][1])
        pc = self._ps[grid_number].copy()

        corners = np.zeros((4,3))

        corners[0,:] = pc
        corners[1,:] = pc + s_side
        corners[2,:] = pc + f_side
        corners[3,:] = pc + s_side + f_side

        return corners


    @property
    def xyz(self):
        return self.to_explicit()


    def to_explicit(self):
        """
        Return the entire grid as an n x 3 array, defining the x,y,z positions
        of each pixel.

        Returns
        -------
        xyz : np.ndarray, float
            An N x 3 array of the x,y,z positions of each pixel. Note that this
            is a flattened version of what you get for each grid individually
            using `grid_as_explicit`.

        See Also
        --------
        grid_as_explicit
        """
        ex_grids = [ self.grid_as_explicit(i) for i in range(self.num_grids) ]
        xyz = np.concatenate([ g.reshape((g.shape[0]* g.shape[1], 3)) for g in ex_grids ])
        return xyz


    def grid_as_explicit(self, grid_number):
        """
        Get the x,y,z coordiantes for a single grid.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get.

        Returns
        -------
        xyz : np.ndarray, float
            An (shape) x 3 array of the x,y,z positions of each pixel

        See Also
        --------
        to_explicit
        """

        p, s, f, shape = self.get_grid(grid_number)

        # xyz = i * s + j * f, where i,j are ints running over range `shape`
        mg = np.mgrid[0:shape[0]-1:1j*shape[0], 0:shape[1]-1:1j*shape[1]]
        xyz = np.outer(mg[0].flatten(), s) + np.outer(mg[1].flatten(), f)
        xyz += p # translate
        xyz = xyz.reshape( (shape[0], shape[1], 3) )

        return xyz

    
    def as_array(self):
        """
        Return the entire BasisGrid object as a 2D numpy array.
        
        Returns
        -------
        psfs : np.ndarray, float
            A 2D array where the first axis is the grid number, the second
            is p_x/p_y/p_z/s_x/s_y/s_z/f_x/f_y/f_z/shape_s/shape_f
        """
        
        psfs = np.zeros((self.num_grids, 11))
        
        for g in range(self.num_grids):
            psfs[g] = np.concatenate(self.get_grid(g))
        
        return psfs
        
        
    @classmethod
    def from_array(cls, arr):
        """
        Create a BasisGrid object from an appropriate 2D numpy array.
        
        Parameters
        ----------
        psfs : np.ndarray, float
            A 2D array where the first axis is the grid number, the second
            is p_x/p_y/p_z/s_x/s_y/s_z/f_x/f_y/f_z/shape_s/shape_f
        """
        
        if not arr.shape[1] == 11:
            raise ValueError('`arr` argument must be shape N x 11')
        
        instance = cls()
        for g in range(arr.shape[0]):
            instance.add_grid(arr[g,0:3],
                              arr[g,3:6],
                              arr[g,6:9],
                              arr[g,9:11].astype(np.int))
        
        return instance