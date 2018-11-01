

"""
Conversion utilities for CSPAD gainmaps.

There are no set conventions for gain maps.... UNFORTUNATELY

WE ADOPT THE FOLLOWING CONVENTIONs:
    -- the gain map for a cspad will be represented by a (32,185,388) array
       with the same readout order as the intensity data (makes sense...)
    -- the values of the elements of this array give the GAIN, so:

          data / gain_map = intensity

You can generate gainmaps easily using the geometry utilities in psgeom, then
save them in DAQ or cheetah format using these utilities.

Note the different conventions known:

-- psgeom / psana
    shape:     (32, 185, 388)
    values:    are the gain value (usually 1.0 for low-, 7.2 for high-gain)

-- cheetah
    shape:     (1480, 1552) [cheetah format]
    values:    are the gain value (usually 1.0 for low-, 7.2 for high-gain)

-- DAQ
    shape:     (11840, 194) [stacked ASICs]
    values:    0 is highgain, 1 is lowgain
"""


import h5py
import numpy as np

from psgeom import translate


def load_cheetah(filename):
    """
    Load a DAQ-formatted gain map.
    
    Parameters
    ----------
    filename : str
        The location on disk to load the gainmap
    
    Returns
    -------
    gainmap : np.ndarray
        psana-formatted gainmap, which is an array of shape (32, 185, 388),
        float value is the gain
    """
    f = h5py.File(filename, 'r')
    gainmap = np.array(f['data/data'])
    f.close()
    return translate_cheetah(gainmap)
    

def translate_cheetah(gainmap):
    """
    Translate a gain map from Cheetah format to DAQ format
    
    Parameters
    ----------
    gainmap : np.ndarray
        Cheeath-formatted gainmap, which is an array of floats of shape
        (1480, 1552). The float value is the gain (usually 1 / ~7 for a)
        CSPAD.
    
    Returns
    -------
    new_gainmap : np.ndarray
        psana-formatted gainmap, which is an array of shape (32, 185, 388),
        float value is the gain
    """
    
    if gainmap.shape != (1480, 1552):
        raise ValueError('gain map is not the expected cheetah shape. '
                         'Expect (1480, 1552), got: %s' % str(gainmap.shape))

    new_gainmap = translate._cheetah_to_twobyones(gainmap)
    
    return new_gainmap
    
    
def write_cheetah(filename, gainmap):
    """
    Write a cheetah-formatted gain map to disk.
    
    Parameters
    ----------
    filename : str
        The location on disk to write the gainmap
    
    gainmap : np.ndarray
        A shape-(32,185,388) gainmap [psana format]
    """
    
    if not gainmap.shape == (32, 185, 388):
        raise ValueError('`gainmap` has wrong shape, expected (32, 185, 388),'
                         ' got: %s' % str(gainmap.shape))
    
    cheetah_image = np.zeros((1480, 1552), dtype=gainmap.dtype)
    
    for q in range(4):
        for a in range(8):

            x_start = 388 * q
            x_stop  = 388 * (q+1)

            y_start = 185 * a
            y_stop  = 185 * (a + 1)

            ng_idx = q * 8 + a
            cheetah_image[y_start:y_stop,x_start:x_stop] = gainmap[ng_idx,:,:]
            
    f = h5py.File(filename, 'w')
    f['/data/data'] = cheetah_image
    f.close()
    print(('Wrote cheetah formatted gainmap: %s' % filename))
    
    return
    

def load_daq(filename, gain_ratio=7.2):
    """
    Load a DAQ-formatted gain map.
    
    Parameters
    ----------
    filename : str
        The location on disk to load the gainmap
    
    gain_ratio : float
        The gain ratio between high- and low-gain mode. IE, ratio of the 
        number of ADU recorded per photon in high gain mode over the same
        value in low gain mode.
    
    Returns
    -------
    gainmap : np.ndarray
        psana-formatted gainmap, which is an array of shape (32, 185, 388),
        float value is the gain
    """
    gainmap = np.loadtxt(filename)
    return translate_daq(gainmap, gain_ratio=gain_ratio)


def translate_daq(gainmap, gain_ratio=7.2):
    """
    Translate a DAQ-formatted gain map to psana format.
    
    Parameters
    ----------
    gainmap : np.ndarray
        DAQ-formatted gainmap, which is an array of shape (11840, 194)
        representing a stack of ASICs. In the DAQ format, high gain
        is 0.0 and low gain is 1.0 (dont ask me...)
    
    gain_ratio : float
        The gain ratio between high- and low-gain mode. IE, ratio of the 
        number of ADU recorded per photon in high gain mode over the same
        value in low gain mode.
    
    Returns
    -------
    new_gainmap : np.ndarray
        psana-formatted gainmap, which is an array of shape (32, 185, 388),
        float value is the gain
    """
    
    if gainmap.shape != (11840, 194):
        raise ValueError('gain map is not the expected DAQ shape. '
                         'Expect (11840, 194), got: %s' % str(gainmap.shape))
        
    new_gainmap = np.zeros((32, 185, 388), dtype=np.float32)
    
    for q in range(4):
        for a in range(8): # which 2x1

            x_start = 388 * q
            x_stop  = 388 * (q + 1)

            y_start = 185 * a
            y_stop  = 185 * (a + 1)

            asic_index = q * 16 + a * 2

            two_by_one = np.concatenate([gainmap[ asic_index   *185:(asic_index+1)*185,:], 
                                         gainmap[(asic_index+1)*185:(asic_index+2)*185,:]], axis=1)
            
            ng_idx = q * 8 + a
            new_gainmap[ng_idx,:,:] = two_by_one
      
    # in DAQ gain maps, 0 is highgain, 1 is lowgain
    highgain = (np.abs(new_gainmap) < 1e-8)
    new_gainmap[highgain] = gain_ratio
    new_gainmap[np.logical_not(highgain)] = 1.0
            
    return new_gainmap


def write_daq(filename, gainmap):
    """
    Write a cheetah-formatted gain map to disk.
    
    Parameters
    ----------
    filename : str
        The location on disk to write the gainmap
    
    gainmap : np.ndarray
        A shape-(32,185,388) gainmap [psana format]
    """
    
    if not gainmap.shape == (32, 185, 388):
        raise ValueError('`gainmap` has wrong shape, expected (32, 185, 388),'
                         ' got: %s' % str(gainmap.shape))
    
    new = np.zeros((11840, 194))

    for i in range(4):
        for j in range(8):
            
            psind = i*8 + j
            sec1, sec2 = np.hsplit(gainmap[psind,:,:], 2)
            
            block = i*16 + j*2
            new[block*185:(block+1)*185]     = sec1
            new[(block+1)*185:(block+2)*185] = sec2
            
    # ensure we get the correct gain values
    # in DAQ format, highgain is 0, lowgain is 1
    highgain = (new > 1.0)
    new[highgain] = 0.0
    new[np.logical_not(highgain)] = 1.0
    
    np.savetxt(filename, new)
    print(('Wrote DAQ formatted gainmap: %s' % filename))
    
    return
    
    
    
