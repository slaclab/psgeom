
"""
Code for fitting data and generating geometries from them. E.g. if you have
geometries at 3 z-motor positions and want to determine some intermediate value.
"""

import numpy as np
from psgeom import basisgrid


class BasisGridInterpolator(object):
    
    def __init__(self, basisgrids, independent_variables):
        """

        Parameters
        ----------
        basisgrids : list of 


        independent_variables : np.ndarray
            Array of k (rows) independent variables measured at n (columns) 
            positions that modify the geometry some how. These are what will
            be interpolated against. The cannonical example would be the
            z-motor position that drives the camera forward and backward.
        """
        
        self.basisgrids = basisgrids
        
        if len(independent_variables.shape) == 1:
            independent_variables = independent_variables.reshape(-1, 1)
        self.independent_variables = independent_variables
        
        # make sure num elements is the same for each grid
        for bg in self.basisgrids:
            if not bg.num_grids == self.basisgrids[0].num_grids:
                raise ValueError('basisgrids must have same number of grids!')
                
        # make sure the shapes of each grid element in for each camera is same
        for g in range(self.grids_per_basis):
            for bg in self.basisgrids:
                if not (bg.get_grid(g)[3] == self.basisgrids[0].get_grid(g)[3]):
                    raise ValueError('grids must be same shape')
        
        
        self._coefficient_matrix = self._interpolate_basis_grids()
        
        return
        
        
    def predict(self, new_position):
        """
        Predict a new geometry via interpolation (or, if ballsy, exterpolation).
        
        Parameters
        ----------
        new_position : np.ndarray
            A vector of positions. Length should match self.num_indept_vars.
            
        Returns
        -------
        bg : basisgrid.BasisGrid
            A basisgrid object describing the estimated geometry
        """
        
        if type(new_position) in [int, float]:
            new_position = np.array(new_position)
        
        if not new_position.shape == (self.num_indept_vars,):
            raise ValueError('`new_position`s shape does not match original '
                             '`independent_variables` input')
                             
                             
        # augment the passed value with a one for the const offsets
        anp = np.hstack([ new_position, np.ones(new_position.shape[0]) ])
        
        # A X --> Y, forward direction
        bg_vector = np.dot(anp, self._coefficient_matrix)
        
        # get bg_vector into bg format
        shps = self.basisgrids[0].as_array()[:,9:] # grid shapes same for all
        psf = bg_vector.reshape(self.grids_per_basis, -1)
        assert psf.shape[0] == shps.shape[0]
        
        bg = basisgrid.BasisGrid.from_array( np.hstack([psf, shps]) )
                             
        return bg
        

    @property
    def num_basis(self):
        return len(self.basisgrids)
    

    @property
    def grids_per_basis(self):
        return self.basisgrids[0].num_grids
        
       
    @property 
    def num_indept_vars(self):
        return self.independent_variables.shape[1]


    @property
    def p_slopes(self):
        """
        Return the infered slopes for each p-vector, final shape is (parameter, xyz)
        """
        px = self._coefficient_matrix[:,0::9]
        py = self._coefficient_matrix[:,1::9]
        pz = self._coefficient_matrix[:,2::9]
        return np.vstack([px, py, pz]).T
    

    def _interpolate_basis_grids(self):
        """
        Perform the LSQ optimization. The final system we want is:
        
                A X = Y
            
            A : n   x k+1,  motor positions (and a constant col) for each postn.
            X : k+1 x 9gn,  coefficient matrix 
            Y : n   x 9gn,  basis grid elements
            
            n : number of unique positions/geometries
            g : number of basis grids
        
        """
        
        # motor positions, and a column of 1's for a constant
        ivs = self.independent_variables.shape
        A = np.hstack([ self.independent_variables, 
                        np.ones((ivs[0], 1)) ])
        
        # measured basis grid elements corresponding to those positions
        Y = np.vstack([ grid.as_array()[:,:9].flatten() \
                        for grid in self.basisgrids ])
        
        assert A.shape == (ivs[0], ivs[1]+1), 'ivs: %s' % (str(ivs),)
        assert Y.shape == (ivs[0], 9 * self.grids_per_basis)
        
        X, resid, rank, s = np.linalg.lstsq(A, Y)
        assert X.shape == (ivs[1]+1, 9 * self.grids_per_basis)
        
        if np.sum(np.abs(resid)) > 1e-8:
            print('Warning: high fitting residuals')
    
        return X
    
    
