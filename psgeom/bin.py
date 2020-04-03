
import numpy as np


class Averager(object):
    
    def __init__(self, point_values, mask, n_bins=101):
        """
        Parameters
        ----------
        point_values : np.ndarray (float)
            For each pixel, this is the point value of that pixel that determines what
            bin it gets assigned to
        mask : np.ndarray (int)
            A boolean (int) saying if each pixel is masked or not
        n_bins : int
            The number of bins to employ. If `None` guesses a good value.
        """

        if not point_values.shape == mask.shape:
            raise ValueError('`point_values` and `mask` must have same shape')
        self.point_values = point_values
        self.mask = mask
        
        self.n_bins = n_bins
        
        # figure out the number of bins to use
        if n_bins != None:
            self.n_bins = n_bins
            self._bin_factor = float(self.n_bins-0.5) / self.point_values.max()
        else:
            self._bin_factor = 25.0
            self.n_bins = (self.point_values.max() * self._bin_factor) + 1
        
        self._bin_assignments = np.floor( point_values * self._bin_factor ).astype(np.int32)
        self._normalization_array = (np.bincount( self._bin_assignments.flatten(),
                                     weights=self.mask.flatten() ) \
                                     + 1e-100).astype(np.float)

        assert self.n_bins == self._bin_assignments.max() + 1, 'bin mismatch in init'
        self._normalization_array = self._normalization_array[:self.n_bins]
        
        return

    
    def __call__(self, image):
        """
        Bin pixels by their point_values
        
        Parameters
        ----------            
        image : np.ndarray
            The intensity at each pixel, same shape as pixel_pos


        Returns
        -------
        bin_centers : ndarray, float
            The center of each bin.

        bin_values : ndarray, int
            The average in the bin.
        """

        if not (image.shape == self.point_values.shape):
            raise ValueError('`image` and `point_values` must have the same shape')

        weights = image.flatten() * self.mask.flatten()
        bin_values = np.bincount(self._bin_assignments.flatten(), weights=weights)
        bin_values /= self._normalization_array
   
        assert bin_values.shape[0] == self.n_bins, 'bin number mismatch (%d, %d)' \
                                                   % (bin_values.shape[0], self.n_bins)
    
        return bin_values

    @property
    def bin_centers(self):
        return np.arange(self.n_bins) / self._bin_factor


