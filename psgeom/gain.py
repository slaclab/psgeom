

def translate_gain_cheetah_to_DAQ(gainmap):
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
        DAQ-formatted gainmap, which is an array of shape (11840, 194)
        representing a stack of ASICs. In the DAQ format, high gain
        is 0.0 and low gain is 1.0 (dont ask me...)
    """
    
    if gainmap.shape != (1480, 1552):
        raise ValueError('gain map is not the expected cheetah shape. '
                         'Expect (1480, 1552), got: %s' % str(gainmap.shape))

    asic_format = translate._cheetah_to_asics(gainmap)
    new = np.zeros((11840, 194))

    for i in range(4):
        for j in range(16):
            block = i * 16 + j
            new[block*185:(block+1)*185] = asic_format[i,j,:,:]

    # in DAQ gain maps, 0 is highgain, 1 is lowgain
    highgain = (new > 1.0)
    new[highgain] = 0.0
    new[np.logical_not(highgain)] = 1.0
    
    return new

def translate_gain_DAQ_to_cheetah(gainmap, gain_ratio=7.2):
    """
    Translate a gain map from Cheetah format to DAQ format
    
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
        Cheeath-formatted gainmap, which is an array of floats of shape
        (1480, 1552). The float value is the gain (usually 1 / ~7 for a)
        CSPAD.
    """
    
    if gainmap.shape != (11840, 194):
        raise ValueError('gain map is not the expected cheetah shape. '
                         'Expect (11840, 194), got: %s' % str(gainmap.shape))
        
    cheetah_map = np.zeros((1480, 1552), dtype=np.float32)
    
    for q in range(4):
        for a in range(8): # which 2x1

            x_start = 388 * q
            x_stop  = 388 * (q + 1)

            y_start = 185 * a
            y_stop  = 185 * (a + 1)

            asic_index = q * 16 + a * 2

            two_by_one = np.concatenate([gainmap[ asic_index   *185:(asic_index+1)*185,:], 
                                         gainmap[(asic_index+1)*185:(asic_index+2)*185,:]], axis=1)
            cheetah_map[y_start:y_stop,x_start:x_stop] = two_by_one
      
    # in DAQ gain maps, 0 is highgain, 1 is lowgain
    highgain = (cheetah_map == 0.0)
    cheetah_map[highgain] = 7.1999998
    cheetah_map[np.logical_not(highgain)] = 1.0
            
    return cheetah_map
