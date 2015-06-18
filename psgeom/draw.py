

"""
A utility library for visualizing geometries
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches

quad_colors = ['k', 'g', 'purple', 'b']


def sketch_2x1s(pixel_positions, mpl_axes=None):
    """
    Draw a rough sketch of the layout of the CSPAD
    Parameters
    ----------
    pixel_positions : np.ndarray
        The x,y,z coordinates of the pixels on the CSPAD
    """
    
    if pixel_positions.shape not in [(4,8,185,388,3), (4,8,185,388,3)]:
        raise ValueError('`pixel_positions` has incorrect shape: '
                         '%s' % str(pixel_positions.shape))

    if not mpl_axes:
        plt.figure()
        ax = plt.subplot(111)
    else:
        ax = mpl_axes

    for i in range(4):
        for j in range(8):
            x = pixel_positions[i,j,:,:,0]
            y = pixel_positions[i,j,:,:,1]
            corners = np.zeros((5,2))

            corners[0,:] = np.array([ x[0,0],   y[0,0] ])     # bottom left
            corners[1,:] = np.array([ x[0,-1],  y[0,-1] ])    # bottom right
            corners[3,:] = np.array([ x[-1,0],  y[-1,0] ])    # top left
            corners[2,:] = np.array([ x[-1,-1], y[-1,-1] ])   # top right
            corners[4,:] = np.array([ x[0,0],   y[0,0] ])     # make rectangle

            ax.plot(corners[:,0], corners[:,1], lw=2, color=quad_colors[i])
            ax.scatter(x[0,0], y[0,0])
            
    beam_center = plt_patches.Circle((0, 0), 2, fill=True, lw=1, color='orange')
    ax.add_patch(beam_center)
            
    # mirror x axis for CXI convention
    if not ax.xaxis_inverted():
        ax.invert_xaxis()

    if mpl_axes:
        return ax
    else:
        plt.show()
