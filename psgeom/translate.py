

"""
translate.py

Translate between different geometry file formats.

-- psana
extension:  .data
format:     flat text
units:      microns


-- Cheetah
extension:   .h5
format:      HDF5
units:       mm for x/y, microns for z


-- CrystFEL
extension:   .geom
format:      flat text
units:       all CrystFEL units are pixel units, 
             except the sample-to-detector offset
references:  http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html


-- Thor
extension:   .dtc
format:      HDF5
units:       intrinsic -- that is, any unit is allowed so long as it is
             self-consistent (all units the same)


-- DIALS
extension:   .json
format       JSON
units:       intrinsic
references:  http://journals.iucr.org/j/issues/2014/04/00/jo5001/index.html
             http://journals.iucr.org/d/issues/2016/04/00/gm5043/index.html#BB30
             http://journals.iucr.org/d/issues/2018/09/00/lp5037/index.html
"""


import os
import re
import getpass
import datetime
import h5py
import math
import warnings
import json

import numpy as np

from psgeom import sensors
from psgeom import basisgrid


def _check_obj(obj):
    """
    check that the object the a detector geometry is being loaded into is "valid"
    """
    # NOT IMPLEMENTED
    return


def _natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


# ---- psana -------------------------------------------------------------------

# define a "type map" that maps a list of known object identifier strings to
# the corresponding types

def map_type(element_name, version=None):
    
    if version is None:
        try:
            version = int( re.search(r'V(\d+)', element_name).group(1) )
        except:
            # most element's version 1 has no "Vx" string
            # will be "MTRX:X1:Y1:X2:Y2" instead of "MTRX:Vx:X1:Y1:X2:Y2"
            version = 1

    if element_name.startswith('SENS2X1'):
        element_type = sensors.Cspad2x1
        
    elif element_name.startswith('PNCCD'):
        element_type = sensors.PnccdQuad
        if version < 2:
            raise TypeError('psgeom v1.0+ requires PNCCD:V2+, got V%d' % version)
        
    elif element_name.startswith('MTRX'):
        if version >= 2:
            element_type = sensors.Mtrx
        elif version == 1:
            element_type = sensors.MtrxV1
        
    elif element_name.startswith('JUNGFRAU'):
        element_type = sensors.JungfrauSegment

    elif element_name.startswith('EPIX10KA'):
        element_type = sensors.Epix10kaSegment

    return element_type, element_name



def load_psana(obj, filename):
    """
    Load a geometry in psana format.
    
    Parameters
    ----------
    filename : str
        The path of the file on disk.
        
    Returns
    -------
    root : detector.CompoundCamera
        The CompoundCamera instance
        
    References
    ----------
    ..[1] https://confluence.slac.stanford.edu/display/PSDM/Detector+Geometry
    """
    
    _check_obj(obj)
    
    print('Loading: %s' % filename)

    # ---- load information into 2 tables: id_info (names) & trt (data)

    # types/names and id numbers -- unique name is 2 fields
    id_info = np.genfromtxt(filename, dtype=np.str, 
                            usecols=(0,1,2,3), autostrip=True)
    #print(id_info)

    # translation & rotation table (trt)
    trt = np.genfromtxt(filename, dtype=np.float64, usecols=list(range(4,13)), autostrip=True)
    assert trt.shape[1] == 9
    translations = trt[:,0:3]
    rotations    = trt[:,3:6] + trt[:,6:9] # just combine rotations & tilts

    n_rows = len(id_info)
    assert trt.shape[0] == n_rows, 'name/data rows dont match (%d/%d)' % \
                                                    (trt.shape[0], n_rows)


    # ---- find the root of the tree, not to be guarenteed first line
    possible_root_rows = list(range(n_rows)) # elimination

    for i in range(n_rows):
        for j in range(n_rows):
        
            # check to see if row i is listed as a child anywhere
            if (id_info[i][0], id_info[i][1]) == (id_info[j][2], id_info[j][3]):
            
                # if so, we know row i is not a root
                #print "%d is child in row %d" % (i,j)
                if i in possible_root_rows:
                    possible_root_rows.remove(i)
            
        
    # check to make sure we have only a single root
    if len(possible_root_rows) != 1:
        raise IOError('Ambiguous tree structure in geometry file. '
                      '%d roots found.' % len(possible_root_rows))
    else:
        root_index = possible_root_rows[0]
        #if id_info[root_index][0] != 'IP':
        #    warnings.warn('Root object is not labeled "IP".')
    

    # ---- traverse tree, adding children / depth frist

    def add_to_tree(current_node_index, parent=None):
    
        cni = current_node_index # shorthand
        cid = (id_info[cni][2], id_info[cni][3])
    
        child_indices = []
        for i in range(n_rows):
            if (id_info[i][0], id_info[i][1]) == cid:
                child_indices.append(i)
        #print '%d --> %s' % (cni, child_indices)
        

        # > if no children, is a SensorElement
        if len(child_indices) == 0:

            # following commented lines deprecated 
            # loop up what type of sensor element we have
            #try:
            #    typ = sensors.type_map[id_info[cni][2]]
            #except KeyError:
            #    raise KeyError('Sensor type: %s not understood.' % id_info[i][2])

            typ, name = map_type(id_info[cni][2])
            
            # TJL note to self:
            # this next line could be problematic if we don't restrict
            # the __init__ method of SensorElements.....
        
            curr = typ.from_type(type_name=name,
                                 id_num=int(id_info[cni][3]),
                                 parent=parent,
                                 rotation_angles=rotations[cni], 
                                 translation=translations[cni])
        
        # > else, is a CompoundCamera
        else:
            curr = obj(type_name=id_info[cni][2],
                       id_num=int(id_info[cni][3]),
                       parent=parent,
                       rotation_angles=rotations[cni], 
                       translation=translations[cni])
                   
            for c in child_indices:
                _ = add_to_tree(c, parent=curr)
            
        # we discard all intermediate results, only root is left
        return curr
            
    root_object = add_to_tree(root_index, parent=None)

    return root_object
    
    
def _mikhail_ordering(list_of_lines):
    """
    for legacy! aka, hopefully we can remove this confusing code soon...
    """
    
    # ordering we're going for starts with all the sensor elements, followed
    # by the quads, followed by whatever else
    
    # assume that the ordering we're dealing with is depth-first
    
    sensors = []
    quads   = []
    other   = []
    root    = []
    
    for line in list_of_lines:
        if 'QUAD' in line[16:]:
            quads.append(line)
        elif 'SENS' in line[16:]:
            sensors.append(line)
        elif 'root' in line[16:]:
            root.append(line)
        else:
            other.append(line)
    
    ordered_list = sensors + quads + other + root
    
    assert len(ordered_list) == len(list_of_lines)
    for e in list_of_lines:
        assert e in ordered_list
    
    return ordered_list
    
    
    
def write_psana(detector, filename, title='geometry'):
    """
    Write a geometry in psana format.

    Parameters
    ----------
    filename : str
        The path of the file on disk.
   
    References
    ----------
    ..[1] https://confluence.slac.stanford.edu/display/PSDM/Detector+Geometry
    """
    
    f = open(filename, 'w')
    
    # write a header
    author = getpass.getuser()
    today = datetime.date.today()
    date = '%d-%d-%d' % (today.month, today.day, today.year)
    time = datetime.datetime.now().time().isoformat()

    header = """
# TITLE       %s
# AUTHOR      %s
# CALIB_TYPE  geometry
# COMMENT:01  WRITTEN BY USER 
# DATE_TIME   %s\t%s
# PARAM:01    PARENT     - name and version of the parent object
# PARAM:02    PARENT_IND - index of the parent object
# PARAM:03    OBJECT     - name and version of the object
# PARAM:04    OBJECT_IND - index of the new object
# PARAM:05    X0         - x-coordinate [um] of the object origin in the parent frame
# PARAM:06    Y0         - y-coordinate [um] of the object origin in the parent frame
# PARAM:07    Z0         - z-coordinate [um] of the object origin in the parent frame
# PARAM:08    ROT_Z      - object design rotation angle [deg] around Z axis of the parent frame
# PARAM:09    ROT_Y      - object design rotation angle [deg] around Y axis of the parent frame
# PARAM:10    ROT_X      - object design rotation angle [deg] around X axis of the parent frame
# PARAM:11    TILT_Z     - object tilt angle [deg] around Z axis of the parent frame
# PARAM:12    TILT_Y     - object tilt angle [deg] around Y axis of the parent frame
# PARAM:13    TILT_X     - object tilt angle [deg] around X axis of the parent frame
# HDR         PARENT IND        OBJECT IND     X0[um]   Y0[um]   Z0[um]   ROT-Z ROT-Y ROT-X     TILT-Z   TILT-Y   TILT-X
""" % (title, author, date, time)

    f.write(header + '\n')

    lines = [] # container for data lines

    fmt_line = '%12s %4d     %12s %4d' + ' %12.6f'*9 + '\n'

    dist = detector.xyz.flatten()[2]

    # write a line for the root node
    root_data = ['IP', 0, detector.type_name, detector.id] + [0.0]*2 + [dist] + [0.0]*6
    root_line = fmt_line % tuple(root_data)
    
    lines.append(root_line)

    # write a line for each child node in the CompoundCamera tree
    def write_children(node):
    
        if hasattr(node, 'children'):
    
            for child in node.children:
                child_data = [node.type_name,
                              node.id,
                              child.type_name,
                              child.id]
                child_data += list(child.translation)
                child_data += list(child.rotation_angles)
                child_data += [0.0]*3
        
                assert len(child_data) == 13
            
                # Set all child z position to 0
                child_data[6] = 0

                line = fmt_line % tuple(child_data)
                lines.append(line)
        
                write_children(child)

    write_children(detector)
   
    # temporary -- for compatability with legacy code -- todo
    # flip the ordering of the lines so that the sensor elements come first,
    # as a lot of existing code requires this ordering
    for l in _mikhail_ordering(lines):
        f.write(l)

    f.close()
    
    return
    

# ---- cheetah -----------------------------------------------------------------
    
    
def _cheetah_to_asics(cheetah_image):
    
    new_image = np.zeros((4,16,185,194), dtype=cheetah_image.dtype)
    
    shape = (185, 388)
    
    num_quads = int(cheetah_image.shape[1] / shape[1])
    if cheetah_image.shape[1] % shape[1] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of quads.' % str(cheetah_image.shape))
    
    num_twoXones = int(cheetah_image.shape[0] / shape[0])
    if cheetah_image.shape[0] % shape[0] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of two-by-ones.' % str(cheetah_image.shape))
    
    for q in range(num_quads):
        for twoXone in range(num_twoXones):
            
            x_start = 388 * q
            x_stop  = 388 * (q+1)
            
            y_start = 185 * twoXone
            y_stop  = 185 * (twoXone + 1)
            
            sec1, sec2 = np.hsplit(cheetah_image[y_start:y_stop,
                                                 x_start:x_stop], 2)
            
            new_image[q,twoXone*2,:,:]   = sec1
            new_image[q,twoXone*2+1,:,:] = sec2
            
            
    return new_image
    
    
def _cheetah_to_twobyones(cheetah_image):
    
    shape = (185, 388)
    
    num_quads = int(cheetah_image.shape[1] / shape[1])
    if cheetah_image.shape[1] % shape[1] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of quads.' % str(cheetah_image.shape))
    
    num_twoXones = int(cheetah_image.shape[0] / shape[0])
    if cheetah_image.shape[0] % shape[0] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of two-by-ones.' % str(cheetah_image.shape))
    
    new_image = np.zeros((num_quads*num_twoXones,185,388), 
                         dtype=cheetah_image.dtype)
    
    for q in range(num_quads):
        for twoXone in range(num_twoXones):
            
            x_start = 388 * q
            x_stop  = 388 * (q+1)
            
            y_start = 185 * twoXone
            y_stop  = 185 * (twoXone + 1)
            
            new_image[q*8 + twoXone,:,:] = cheetah_image[y_start:y_stop,
                                                         x_start:x_stop]
            
    return new_image
    
    
def load_cheetah(obj, filename, pixel_size=109.92, element_type=sensors.Mtrx):
    """
    Load a geometry in cheetah format.
    
    Parameters
    ----------
    filename : str
        The path of the file on disk.
        
    Returns
    -------
    root : camera.Cspad
        The Cspad instance
    """
    
    _check_obj(obj)

    f = h5py.File(filename, 'r')

    if not list(f.keys()) == ['x', 'y', 'z']:
        raise IOError('File: %s is not a valid pixel map, should contain fields'
                      ' ["x", "y", "z"] exlusively' % filename)

    cheetah_shape = f['x'].shape

    # convert m --> um, ends up not mattering tho...
    # also flip the sign of x : cheetah uses +x away from hutch door
    x = -1.0 * _cheetah_to_asics( np.array(f['x']) * 1000000.0 )
    y =        _cheetah_to_asics( np.array(f['y']) * 1000000.0 )

    # for some reason z is in microns, so leave it
    z = _cheetah_to_asics( np.array(f['z']) )

    f.close()

    bg = basisgrid.BasisGrid()
    shape = (185, 388) # will always be this for each asic
    
    # find out how many quads/asics we expect based on the size of the maps
    # we need this logic here to deal with both big and 2x2 CSPADs
    num_quads = int(cheetah_shape[1] / shape[1])
    if cheetah_shape[1] % shape[1] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of quads.' % str(cheetah_shape))
    
    num_asics = int(cheetah_shape[0] / shape[0]) * 2
    if cheetah_shape[0] % shape[0] != 0:
        raise IOError('Unexpected geometry array shape: %s. Could not infer '
                      'number of two-by-ones.' % str(cheetah_shape))

    # loop over each asic, and convert it into a basis grid
    for i in range(num_quads):
        for j in range(num_asics):

            # extract all the corner positions (code ineligant but explicit)
            # corners are numbered 0 -> 4, starting top left and continuing cw
            corners = np.zeros(( 4, 3 ))
            corners[0,:] = ( x[i,j,0,0],   y[i,j,0,0],   z[i,j,0,0]   )
            corners[1,:] = ( x[i,j,0,-1],  y[i,j,0,-1],  z[i,j,0,-1]  )
            corners[2,:] = ( x[i,j,-1,-1], y[i,j,-1,-1], z[i,j,-1,-1] )
            corners[3,:] = ( x[i,j,-1,0],  y[i,j,-1,0],  z[i,j,-1,0]  )

            # average the vectors formed by the corners to find f/s vects
            # the fast scan direction is the last index, s is next
            # f points left -> right, s points bottom -> top
            f = (( corners[1,:] - corners[0,:] ) + ( corners[2,:] - corners[3,:] ))
            s = (( corners[3,:] - corners[0,:] ) + ( corners[2,:] - corners[1,:] ))

            # make them pixel-size magnitude
            f = f * (pixel_size / np.linalg.norm(f))
            s = s * (pixel_size / np.linalg.norm(s))
            
            # p is just location of 1st pixel in memory, which is our 1st corner
            p = corners[0,:]
            bg.add_grid(p, s, f, shape)
            
    geom_instance = obj.from_basisgrid(bg, element_type=element_type)
                 
    return geom_instance


def write_cheetah(detector, filename="pixelmap-cheetah-raw.h5"):
    """
    Write a geometry in cheetah format.

    Parameters
    ----------
    detector : camera.Cspad
    
    filename : str
        The path of the file on disk.
    """

    coordinates = ['x', 'y', 'z']
    
    if not hasattr(detector, 'xyz'):
        raise TypeError('passed `detector` object must have an xyz attr')

    pp = np.squeeze(detector.xyz)
    if pp.shape == (4, 8, 185, 388, 3):
        num_quads = 4
        num_twoXones = 8
    elif pp.shape == (2, 185, 388, 3):
        num_quads = 1
        num_twoXones = 2
        pp = pp.reshape(1, 2, 185, 388, 3) # dummy quad axis
    else:
        raise ValueError('Geometry does not appear to be a CSPAD. xyz '
                         'shape: %s' % str(pp.shape))

    # write an h5
    f = h5py.File(filename, 'w')

    # iterate over x/y/z
    for xyz in range(len(coordinates)):

        cheetah_image = np.zeros((1480, 1552), dtype=np.float32)

        # iterate over each 2x1/quad (note switch)
        for q in range(num_quads):
            for a in range(num_twoXones): # which 2x1

                x_start = 388 * q
                x_stop  = 388 * (q+1)

                y_start = 185 * a
                y_stop  = 185 * (a + 1)

                # if x axis, flip sign
                # unless z axis, convert um --> m ; z-axis in cheetah is in um
                if xyz == 0:
                    unit_factor = - 1.0 / 1000000.0
                elif xyz == 1:
                    unit_factor = 1.0 / 1000000.0
                elif xyz == 2:
                    unit_factor = 1.0
                
                cheetah_image[y_start:y_stop,x_start:x_stop] = unit_factor * pp[q,a,:,:,xyz]


        f['/%s' % coordinates[xyz]] = cheetah_image

    f.close()

    return
    

# ---- crystfel --------------------------------------------------------------- 

def load_crystfel(obj, filename, element_type=sensors.Mtrx, verbose=True):
    """
    Convert a CrystFEL geom file to a Cspad object.
    
    NOTE ON UNITS: all CrystFEL units are pixel units, except the
    sample-to-detector offset, which is typically in meters.
    
    
    Parameters
    ----------
    filename : str
        The path of the file on disk.
        
    Returns
    -------
    root : camera.Cspad
        The Cspad instance
    """
    
    _check_obj(obj)
    
    # NOTE ON UNITS: all CrystFEL units are pixel units, except the
    # sample-to-detector offset
    
    if not filename.endswith('.geom'):
        raise IOError('Can only read flat text files with extension `.geom`.'
                      ' Got: %s' % filename)

    if verbose:
        print("Converting geometry in: %s ..." % filename)
        
    f = open(filename, 'r')
    geom_txt = f.read()
    f.close()


    bg = basisgrid.BasisGrid()


    # measure the absolute detector offset
    re_pz_global = re.search(r'\ncoffset\s+=\s+(\d+.\d+)', geom_txt) 
    if re_pz_global == None:
        print("WARNING: Could not find `coffset` field, defaulting z-offset to 0.0")
        p_z_global = 0.0
    else:
        p_z_global = float(re_pz_global.group(1)) * 1e6 # m --> micron
        if verbose:
            print('Found global z-offset (coffset): %f' % p_z_global)


    # figure out the pixel size
    re_pixel_size = re.search(r'\nres\s+=\s+(\d+.\d+)', geom_txt) 
    if re_pixel_size == None:
        pixel_size = None
    else:
        pixel_size = 1e6 / float(re_pixel_size.group(1)) # m --> micron
        if verbose:
            print('Found pixel size (res) [micron]: %f' % pixel_size)
    
    
    # find out which panels we have to look for
    #   > try to do this in a general way by looking at the unique set of items
    #   > with a corner field    
    panels = re.findall(r'(\S+)/corner', geom_txt)
    panels = _natural_sort(list(set( panels )))
    
    if len(panels) == 0:
        raise IOError('Could not find any panels in file: %s' % filename)
    
    # iterate over each quad / ASIC    
    for panel in panels:

        if verbose:
            print("Reading geometry for: %s" % panel)

        try:
            
            # get pixel size on a per-panel basis
            if pixel_size is None:
                re_pixel_size = re.search(r'%s/res\s+=\s+(\d+.\d+)' % panel, geom_txt) 
                if re_pixel_size == None:
                    raise IOError('could not find required `res` field in file')
                else:
                    pixel_size = 1e6 / float(re_pixel_size.group(1)) # m -> um
                    if verbose:
                        print(('Found pixel size for panel %s (res) [micron]: '
                              '%f' % (panel, pixel_size)))
            

            # match f/s vectors (TODO probably should make this a fxn)
            # fs
            re_fs_x = re.search(r'%s/fs\s+=.*?((.)?\d+.\d+)x' % panel, geom_txt)
            re_fs_y = re.search(r'%s/fs\s+=.*?((.)?\d+.\d+)y' % panel, geom_txt)

            if (re_fs_x is None) and (re_fs_y is None):
                raise IOError('could not find valid ``%s/fs`` field' % panel)

            if re_fs_x is None:
                f_x = 0.0
            else:
                f_x = - float( re_fs_x.group(1) )
                
            if re_fs_y is None:
                f_y = 0.0
            else:
                f_y =   float( re_fs_y.group(1) )
                
            f = np.array([f_x, f_y, 0.0])
            f = f * (pixel_size / np.linalg.norm(f))
            
            # ss
            re_ss_x = re.search(r'%s/ss\s+=.*?((.)?\d+.\d+)x' % panel, geom_txt)
            re_ss_y = re.search(r'%s/ss\s+=.*?((.)?\d+.\d+)y' % panel, geom_txt)

            if (re_ss_x is None) and (re_ss_y is None):
                raise IOError('could not find valid ``%s/ss`` field' % panel)

            if re_ss_x is None:
                s_x = 0.0
            else:
                s_x = - float( re_ss_x.group(1) )
                
            if re_ss_y is None:
                s_y = 0.0
            else:
                s_y =   float( re_ss_y.group(1) )
            
            s = np.array([s_x, s_y, 0.0])
            s = s * (pixel_size / np.linalg.norm(s))
            
            
            # min_fs & min_ss
            re_min_fs = re.search(r'%s/min_fs = (\d+)' % panel, geom_txt)
            re_max_fs = re.search(r'%s/max_fs = (\d+)' % panel, geom_txt)
            
            re_min_ss = re.search(r'%s/min_ss = (\d+)' % panel, geom_txt)
            re_max_ss = re.search(r'%s/max_ss = (\d+)' % panel, geom_txt)
            
            shp = ( np.abs(int(re_max_ss.group(1)) - int(re_min_ss.group(1))) + 1, 
                    np.abs(int(re_max_fs.group(1)) - int(re_min_fs.group(1))) + 1)

            sf_angle = np.degrees( np.arccos( np.dot(s, f) / np.square(pixel_size) ) )

            print(panel, sf_angle)
            
        except AttributeError as e:
            print(e)
            raise IOError('Geometry file incomplete -- cant parse one or '
                          'more basis vector fields (ss/fs) for panel: %s' % panel)

        # match corner postions, that become the p vector
        # note we have to convert from pixel units to mm
        # and also that CrystFEL measures the corner from the actual
        # *corner*, and not the center of the corner pixel!
        
        # also, remember the s[0] and f[0] have already been x-flipped
        
        try:
            
            re_cx = re.search(r'%s/corner_x\s+=\s+((.)?\d+(.\d+)?)' % panel, geom_txt)
            p_x = - float( re_cx.group(1) ) * pixel_size + 0.5 * (s[0] + f[0])

            re_cy = re.search(r'%s/corner_y\s+=\s+((.)?\d+(.\d+)?)' % panel, geom_txt)
            p_y =   float( re_cy.group(1) ) * pixel_size + 0.5 * (s[1] + f[1])
            
            
            # it's allowed to also have individual z-offsets for
            # each panel, so look for those (CrystFEL units: meters)                
            re_cz = re.search(r'%s/coffset\s+=\s+((.)?\d+.\d+)' % panel, geom_txt)
            if re_cz == None:
                if verbose:
                    print('Could not find z data for %s' % panel)
                p_z = p_z_global 
            else:
                # add to the global offset
                p_z = p_z_global + float( re_cz.group(1) ) * 1e6 # m --> micron

            p = np.array([p_x, p_y, p_z])

        except AttributeError as e:
            print(e)
            raise IOError('Geometry file incomplete -- cant parse one or '
                          'more corner fields for panel: %s' % panel)

        # finally, add the ASIC to the basis grid
        bg.add_grid(p, s, f, shp)

    if verbose:
        print(" ... successfully converted geometry.")
    
    geom_instance = obj.from_basisgrid(bg, element_type=element_type)
    
    return geom_instance
    

def write_generic_crystfel(detector, filename, coffset=None, **kwargs):
    """
    Parameters
    ----------
    detector : camera.CompoundAreaCamera
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Should end in '.geom'

    coffset: float
        Detector home position to sample distance in metres. 
        When coffset is None, coffset is set to detector distance.
    """
    
    bg = detector.to_basisgrid()
    
    def get_sign(v):
        if v >= 0:
            s = '+'
        else:
            s = '-'
        return s

    with open(filename, 'w') as of:
    
        of.write("; This file contains a geometry generated by psgeom\n")
        of.write("; https://github.com/slaclab/psgeom\n")
    
        of.write(generic_header + '\n')
        
        if 'maskfile' in kwargs: 
            of.write('mask_file = ' + str(kwargs['maskfile']) + '\n')
            of.write('mask = /entry_1/data_1/mask\n')
            of.write('mask_good = 0x0000\n')
            of.write('mask_bad = 0xffff\n')
        else:
            of.write('; mask = /entry_1/data_1/mask\n')
            of.write('; mask_good = 0x0000\n')
            of.write('; mask_bad = 0xffff\n')
        of.write('\n')

        # remember: cheetah data is flattened along the slowest dimension
        # so that (n_elements, slow, fast) --> (n_elements * slow, fast)
        # this final 2d array is called "the slab"

        grid_index  = 0 # count the basisgrid we are on
        ss_slab_pos = 0 # the slow index of "the slab"

        for leaf_index, leaf in enumerate(detector.leaves):

            sp_shape = leaf.subpanel_shape
            for sp_ss_index in range(sp_shape[0]):
                for sp_fs_index in range(sp_shape[1]):

                    subpanel_index = sp_ss_index * sp_shape[1] + sp_fs_index

                    p, s, f, sp = bg.get_grid(grid_index)
                    panel_name = "p%da%d" % (leaf_index, subpanel_index)
                    
                    # write the basis vectors           
                    f_sqt = math.sqrt(f[0]**2 + f[1]**2)
                    s_sqt = math.sqrt(s[0]**2 + s[1]**2)

                    pixel_size = f_sqt
                     
                    of.write("%s/fs = %s%fx %s%fy\n" % ( panel_name,
                                                         get_sign(-f[0]/f_sqt), abs(f[0]/f_sqt), 
                                                         get_sign( f[1]/f_sqt), abs(f[1]/f_sqt) ))
                    of.write("%s/ss = %s%fx %s%fy\n" % ( panel_name,
                                                         get_sign(-s[0]/s_sqt), abs(s[0]/s_sqt), 
                                                         get_sign( s[1]/s_sqt), abs(s[1]/s_sqt) ))
                    of.write("%s/res = %.3f\n" % (panel_name, 1e6 / pixel_size)) # um --> m
                    
                    # write the corner positions
                    tagcx = "%s/corner_x" % panel_name
                    tagcy = "%s/corner_y" % panel_name
                    tagcz = "%s/coffset"  % panel_name
                
                    # CrystFEL measures the corner from the actual *corner*, and not
                    # the center of the corner pixel (dont forget to x-flip s[0], f[0])
                    
                    cx = - float(p[0])/pixel_size + 0.5 * (f[0] + s[0])/pixel_size
                    cy =   float(p[1])/pixel_size - 0.5 * (f[1] + s[1])/pixel_size
                    
                    of.write("%s = %f\n" % (tagcx, cx))
                    of.write("%s = %f\n" % (tagcy, cy))

                    # the z-axis is in *** meters ***
                    if coffset is None:
                        dist = float(p[2]) / 1e6
                    else:
                        dist = coffset
                    of.write("%s = %f\n" % (tagcz, dist))

                    # combine leaf & slow axes
                    #    note  : sp here is just the shape of the basis grid
                    #    note  : sp_ss_index/sp_fs_index is where we are in the subpanel
                    #    note  : fs/ss indices are INCLUSIVE

                    of.write("%s/min_fs = %d\n" % (panel_name, sp[1] * sp_fs_index ))
                    of.write("%s/max_fs = %d\n" % (panel_name, sp[1] * (1 + sp_fs_index) - 1))

                    # ss_slab_pos is incremented below, outside this loop
                    of.write("%s/min_ss = %d\n" % (panel_name, ss_slab_pos))
                    of.write("%s/max_ss = %d\n" % (panel_name, ss_slab_pos + sp[0] - 1))
                    
                    # this tells CrystFEL to use this panel
                    of.write("%s/no_index = 0\n" % panel_name)
                    
                    of.write("\n") # new line

                    grid_index += 1

                # // end fs for loop
                ss_slab_pos += sp[0]

            # // end ss for loop
        # // end leaves for loop
    # // end with open(...) statement

    # make sure things add up :)
    assert grid_index == bg.num_grids, grid_index
            
    return


def write_cspad_crystfel(detector, filename, coffset=None, intensity_file_type='cheetah',
                         pixel_size=109.92, **kwargs):
    """
    Write a CSPAD geometry to disk in CrystFEL format. Note that some fields
    will be written but left blank -- these are fields you probably should
    fill in before performing any computations in CrystFEL, but are information
    that psgeom has no handle on (e.g. detector gain).
    
    Thanks to Rick Kirian & Tom White for assistance with this function.
    
    Parameters
    ----------
    detector : cspad.CSPad
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Should end in '.geom'
    
    coffset: float
        Detector home position to sample distance in metres.
        When coffset is None, coffset is set to detector distance.
    
    Optional Parameters
    -------------------
    intensity_file_type : str, {'cheetah'}
        The kind of file this geometry file will be used with. Necessary to tell
        CrystFEL how intensity data map onto the detector

    pixel_size : float
        Pixel size in microns

    maskfile : str
        Hdf5 filename of a mask used to indexing and integration by CrystFEL.
    """
    
    bg = detector.to_basisgrid()
    
    def get_sign(v):
        if v >= 0:
            s = '+'
        else:
            s = '-'
        return s
    
    
    if intensity_file_type == 'cheetah':
        
        # this is complex, so I went the lazy route and copied an
        # existing file
        intensity_map = crystfel_cheetah_intensities.split('-')
        assert len(intensity_map) == 64
        
    else:
        raise ValueError('Cannot write geometries for '
                        '`intensity_file_type`: %s, only currently '
                        'have implemented writers for '
                        '{"cheetah"}' % intensity_file_type)
    
    
    
    with open(filename, 'w') as of:
    
        of.write("; This file contains a CSPAD geometry generated by psgeom\n")
        of.write("; https://github.com/slaclab/psgeom\n")

        if coffset is None:
            of.write(cspad_header_noClen + '\n')
        else:
            of.write(cspad_header + '\n')

        if 'maskfile' in kwargs: 
            of.write('mask_file = ' + str(kwargs['maskfile']) + '\n')
            of.write('mask = /entry_1/data_1/mask\n')
            of.write('mask_good = 0x0000\n')
            of.write('mask_bad = 0xffff\n')
        else:
            of.write('; mask = /entry_1/data_1/mask\n')
            of.write('; mask_good = 0x0000\n')
            of.write('; mask_bad = 0xffff\n')

        of.write(cspad_groups + '\n')
    
        # iterate over each basis grid object
        # for a full CSPAD, this will be 64 elements
        # for a 2x2, it will be 4 elements and the "quad" will always be 0
        for grid_index in range(bg.num_grids):
                
            asic = grid_index % 16
            quad = grid_index / 16
            
            p, s, f, sp = bg.get_grid(grid_index)

            panel_name = "q%da%d" % (quad, asic)
            
            # tell crystFEL how read intensity values in a file
            of.write(intensity_map[grid_index].strip()  + '\n')
            
            # write the basis vectors           
            sqt = math.sqrt(f[0]**2 + f[1]**2) 
            of.write("%s/fs = %s%fx %s%fy\n" % ( panel_name,
                                                   get_sign(-f[0]/sqt), abs(f[0]/sqt), 
                                                   get_sign( f[1]/sqt), abs(f[1]/sqt) ))
            sqt = math.sqrt(s[0]**2 + s[1]**2)
            of.write("%s/ss = %s%fx %s%fy\n" % ( panel_name,
                                                   get_sign(-s[0]/sqt), abs(s[0]/sqt), 
                                                   get_sign( s[1]/sqt), abs(s[1]/sqt) ))
            
            # write the corner positions
            tagcx = "%s/corner_x" % panel_name
            tagcy = "%s/corner_y" % panel_name
            tagcz = "%s/coffset"  % panel_name
        
            # CrystFEL measures the corner from the actual *corner*, and not
            # the center of the corner pixel (dont forget to x-flip s[0], f[0])
            
            cx = - float(p[0])/pixel_size + 0.5 * (f[0] + s[0])/pixel_size
            cy =   float(p[1])/pixel_size - 0.5 * (f[1] + s[1])/pixel_size
            
            of.write("%s = %f\n" % (tagcx, cx))
            of.write("%s = %f\n" % (tagcy, cy))
            
            # the z-axis is in *** meters ***
            if coffset is None:
                dist = float(p[2]) / 1e6
            else:
                dist = coffset
            of.write("%s = %f\n" % (tagcz, dist ))
            
            # this tells CrystFEL to use this panel
            of.write("%s/no_index = 0\n" % panel_name)
            
            of.write("\n") # new line
    
    return    

# ---- DIALS -------------------------------------------------------------------


def load_dials(obj, filename, scale_factor=1000.0):
    """
    Load a geometry in DIALS format.
    
    Parameters
    ----------
    obj : camera.CompoundAreaCamera
        An instance of CompoundAreaCamera or a subclass of it

    filename : str
        The path of the file on disk.

    scale_factor : float
        A factor to scale the units by (to e.g. meet a new convention for units)
        
    Returns
    -------
    obj : variable
        A populated instance of the camera object
    """

    with open(filename, 'r') as f:
        jo = json.load(f)

    d = jo["detector"]
    if len(d) != 1:
        raise NotImplementedError('multiple detectors in json')
    else:
        base = jo["detector"][0]

    to_visit = []
    bg_tmp = {}

    def get_F(node):
   
        # see Acta Cryst. D (2018). D74, 877-894 eqs (1) & (2)

        dx = np.array(node["fast_axis"]) # x = fast
        dy = np.array(node["slow_axis"]) # y = slow
        dn = np.cross(dx, dy)            # dn = dx x dy
        d0 = np.array(node["origin"])    
        F  = np.array([dx, dy, dn, d0]).T
        F  = np.vstack([F, np.array([0, 0, 0, 1])])

        return F

    def dfs(current_node, cF):

        # the geometry file is separated into two parts, first a
        # 'hierarchy' that specifies the translations/rotations
        # and points to a 'panel' index
        #
        # the second part is the 'panels', which are the leaves
        # of the hierarchy, but are located in a different part
        # of the JSON file -- these have one final trans/rot to
        # apply

        if "children" in list(current_node.keys()):
            for i in range(len(current_node["children"])):
                new_node = (current_node["children"][i],
                            np.dot(cF, get_F(current_node)) )
                to_visit.append( new_node )

        if "panel" in list(current_node.keys()):
            panel_index = current_node["panel"]
            panel = base["panels"][panel_index]
            #print "adding panel:", panel_index

            px_size = np.array(panel["pixel_size"])
            shp = np.array(panel["image_size"])

            pF = get_F(panel)
            final_F = np.dot(cF, pF)

            # extract the final (transformed) panel vectors
            p = final_F[:3,3] * scale_factor
            s = final_F[:3,1] * px_size[1] * scale_factor
            f = final_F[:3,0] * px_size[0] * scale_factor

            # DIAS convention is +z points from IP to source
            # this is the x-flipped version of our convention
            p[0] *= -1
            s[0] *= -1
            f[0] *= -1

            #bg.add_grid(p, s, f, shp)
            bg_tmp[panel_index] = (p, s, f, shp)

        if len(to_visit) > 0:
            dfs(*to_visit.pop())

        return

    dfs(base["hierarchy"], get_F(base["hierarchy"]))

    # this code checks to make sure we got all the panels
    # and that they get added in order
    bg = basisgrid.BasisGrid()
    for k in range(max(bg_tmp.keys())+1):
        if k in list(bg_tmp.keys()):
            bg.add_grid(*bg_tmp[k])
        else:
            print('WARNING: panel %d seems to be missing' % k)
    
    geom_instance = obj.from_basisgrid(bg)

    return geom_instance

    
# ---- generic text ------------------------------------------------------------
    
    
def write_psf_text(detector, filename):
    """
    Write a geometry to disk in the following format:
    
        p_x p_y p_z     s_x s_y s_z     f_x f_y f_z
        ...
        
    and include some comments.    
    
    
    Parameters
    ----------
    geometry : cspad.CSPad
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Will end in '.dtc'
    """
    
    # generate a preamble
    preamble = """
# This file contains a CSPAD geometry generated by psgeom: 
# https://github.com/LinacCoherentLightSource/psgeom
#
# The following is a basis grid representation with the following vectors
#
#   p : position vector for an ASIC
#   s : slow-scan pixel vector
#   f : fast-scan pixel vector
#
# all units are mm. Each ASIC is 185 x 194 pixels.
# See the psgeom documentation for more information.


#            p_x        p_y        p_z           s_x        s_y        s_z           f_x        f_y        f_z
"""

    # loop over each grid element and add it to the file
    bg = detector.to_basisgrid()
    body = ""
    
    
    def format(s, total_len=10):
        """
        A little formatting function
        """
        sf = '%.5f' % s
        pad = total_len - len(sf)
        if pad > 0:
            sf = ' ' * pad + sf
        return sf
        
    
    for i in range(bg.num_grids):
        
        # if we're starting a new quad, note that in the file
        if i % 16 == 0:
            body += ('\n# QUAD %d\n' % (i/16))
        
        # add the basis grid
        p, s, f, shp = bg.get_grid(i)
        strp = ' '.join( [ format(x) for x in p ] )
        strs = ' '.join( [ format(x) for x in s ] )
        strf = ' '.join( [ format(x) for x in f ] )
        
        tb = ' ' * 4
        asic = str(i)
        if len(asic) == 1:
            asic = ' ' + asic
        
        body += (asic + tb + strp + tb + strs + tb + strf + '\n')
        
    f = open(filename, 'w')
    f.write(preamble + body)
    f.close()
    
    print("Wrote CSPAD to text at: %s" % filename)
    
    return
    
    

#  ---------- REFERENCE DATA ---------------------------------------------------

generic_header = """
; --- VALUES YOU MAY WANT TO FILL IN MANUALLY ---
; we cannot guarentee these values are what you desire
; however they are filled in with some decent defaults

; clen =  /LCLS/detector_1/EncoderValue
; photon_energy = /LCLS/photon_energy_eV
; adu_per_eV = 0.1

data = /entry_1/data_1/data

dim0 = %
dim1 = ss
dim2 = fs

peak_list = /entry_1/result_1
"""

generic_header_noClen = """
; --- VALUES YOU MAY WANT TO FILL IN MANUALLY ---
; we cannot guarentee these values are what you desire
; however they are filled in with some decent defaults

; clen =  /LCLS/detector_1/EncoderValue
photon_energy = /LCLS/photon_energy_eV
adu_per_eV = 0.1

data = /entry_1/data_1/data

dim0 = %
dim1 = ss
dim2 = fs
"""


cspad_header = """
; --- VALUES YOU MAY WANT TO FILL IN MANUALLY ---
; we cannot guarantee these values are what you desire
; however they are filled in with some decent defaults
; for the large CSPAD detector.

clen =  /LCLS/detector_1/EncoderValue
photon_energy = /LCLS/photon_energy_eV
res = 9097.52
adu_per_eV = 0.00338

data = /entry_1/data_1/data

dim0 = %
dim1 = ss
dim2 = fs

"""

cspad_header_noClen = """
; --- VALUES YOU MAY WANT TO FILL IN MANUALLY ---
; we cannot guarantee these values are what you desire
; however they are filled in with some decent defaults
; for the large CSPAD detector. Note that clen is 
; commented out.

;clen =  /LCLS/detector_1/EncoderValue
photon_energy = /LCLS/photon_energy_eV
res = 9097.52
adu_per_eV = 0.00338

data = /entry_1/data_1/data

dim0 = %
dim1 = ss
dim2 = fs

"""

cspad_groups = """

; The following lines define "rigid groups" which express the physical
; construction of the detector.  This is used when refining the detector
; geometry.

rigid_group_q0 = q0a0,q0a1,q0a2,q0a3,q0a4,q0a5,q0a6,q0a7,q0a8,q0a9,q0a10,q0a11,q0a12,q0a13,q0a14,q0a15
rigid_group_q1 = q1a0,q1a1,q1a2,q1a3,q1a4,q1a5,q1a6,q1a7,q1a8,q1a9,q1a10,q1a11,q1a12,q1a13,q1a14,q1a15
rigid_group_q2 = q2a0,q2a1,q2a2,q2a3,q2a4,q2a5,q2a6,q2a7,q2a8,q2a9,q2a10,q2a11,q2a12,q2a13,q2a14,q2a15
rigid_group_q3 = q3a0,q3a1,q3a2,q3a3,q3a4,q3a5,q3a6,q3a7,q3a8,q3a9,q3a10,q3a11,q3a12,q3a13,q3a14,q3a15

rigid_group_a0 = q0a0,q0a1
rigid_group_a1 = q0a2,q0a3
rigid_group_a2 = q0a4,q0a5
rigid_group_a3 = q0a6,q0a7
rigid_group_a4 = q0a8,q0a9
rigid_group_a5 = q0a10,q0a11
rigid_group_a6 = q0a12,q0a13
rigid_group_a7 = q0a14,q0a15
rigid_group_a8 = q1a0,q1a1
rigid_group_a9 = q1a2,q1a3
rigid_group_a10 = q1a4,q1a5
rigid_group_a11 = q1a6,q1a7
rigid_group_a12 = q1a8,q1a9
rigid_group_a13 = q1a10,q1a11
rigid_group_a14 = q1a12,q1a13
rigid_group_a15 = q1a14,q1a15
rigid_group_a16 = q2a0,q2a1
rigid_group_a17 = q2a2,q2a3
rigid_group_a18 = q2a4,q2a5
rigid_group_a19 = q2a6,q2a7
rigid_group_a20 = q2a8,q2a9
rigid_group_a21 = q2a10,q2a11
rigid_group_a22 = q2a12,q2a13
rigid_group_a23 = q2a14,q2a15
rigid_group_a24 = q3a0,q3a1
rigid_group_a25 = q3a2,q3a3
rigid_group_a26 = q3a4,q3a5
rigid_group_a27 = q3a6,q3a7
rigid_group_a28 = q3a8,q3a9
rigid_group_a29 = q3a10,q3a11
rigid_group_a30 = q3a12,q3a13
rigid_group_a31 = q3a14,q3a15

rigid_group_collection_quadrants = q0,q1,q2,q3
rigid_group_collection_asics = a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31

; -----------------------------------------------

"""


crystfel_cheetah_intensities = """q0a0/min_fs = 0
q0a0/min_ss = 0
q0a0/max_fs = 193
q0a0/max_ss = 184
-
q0a1/min_fs = 194
q0a1/min_ss = 0
q0a1/max_fs = 387
q0a1/max_ss = 184
-
q0a2/min_fs = 0
q0a2/min_ss = 185
q0a2/max_fs = 193
q0a2/max_ss = 369
-
q0a3/min_fs = 194
q0a3/min_ss = 185
q0a3/max_fs = 387
q0a3/max_ss = 369
-
q0a4/min_fs = 0
q0a4/min_ss = 370
q0a4/max_fs = 193
q0a4/max_ss = 554
-
q0a5/min_fs = 194
q0a5/min_ss = 370
q0a5/max_fs = 387
q0a5/max_ss = 554
-
q0a6/min_fs = 0
q0a6/min_ss = 555
q0a6/max_fs = 193
q0a6/max_ss = 739
-
q0a7/min_fs = 194
q0a7/min_ss = 555
q0a7/max_fs = 387
q0a7/max_ss = 739
-
q0a8/min_fs = 0
q0a8/min_ss = 740
q0a8/max_fs = 193
q0a8/max_ss = 924
-
q0a9/min_fs = 194
q0a9/min_ss = 740
q0a9/max_fs = 387
q0a9/max_ss = 924
-
q0a10/min_fs = 0
q0a10/min_ss = 925
q0a10/max_fs = 193
q0a10/max_ss = 1109
-
q0a11/min_fs = 194
q0a11/min_ss = 925
q0a11/max_fs = 387
q0a11/max_ss = 1109
-
q0a12/min_fs = 0
q0a12/min_ss = 1110
q0a12/max_fs = 193
q0a12/max_ss = 1294
-
q0a13/min_fs = 194
q0a13/min_ss = 1110
q0a13/max_fs = 387
q0a13/max_ss = 1294
-
q0a14/min_fs = 0
q0a14/min_ss = 1295
q0a14/max_fs = 193
q0a14/max_ss = 1479
-
q0a15/min_fs = 194
q0a15/min_ss = 1295
q0a15/max_fs = 387
q0a15/max_ss = 1479
-
q1a0/min_fs = 388
q1a0/min_ss = 0
q1a0/max_fs = 581
q1a0/max_ss = 184
-
q1a1/min_fs = 582
q1a1/min_ss = 0
q1a1/max_fs = 775
q1a1/max_ss = 184
-
q1a2/min_fs = 388
q1a2/min_ss = 185
q1a2/max_fs = 581
q1a2/max_ss = 369
-
q1a3/min_fs = 582
q1a3/min_ss = 185
q1a3/max_fs = 775
q1a3/max_ss = 369
-
q1a4/min_fs = 388
q1a4/min_ss = 370
q1a4/max_fs = 581
q1a4/max_ss = 554
-
q1a5/min_fs = 582
q1a5/min_ss = 370
q1a5/max_fs = 775
q1a5/max_ss = 554
-
q1a6/min_fs = 388
q1a6/min_ss = 555
q1a6/max_fs = 581
q1a6/max_ss = 739
-
q1a7/min_fs = 582
q1a7/min_ss = 555
q1a7/max_fs = 775
q1a7/max_ss = 739
-
q1a8/min_fs = 388
q1a8/min_ss = 740
q1a8/max_fs = 581
q1a8/max_ss = 924
-
q1a9/min_fs = 582
q1a9/min_ss = 740
q1a9/max_fs = 775
q1a9/max_ss = 924
-
q1a10/min_fs = 388
q1a10/min_ss = 925
q1a10/max_fs = 581
q1a10/max_ss = 1109
-
q1a11/min_fs = 582
q1a11/min_ss = 925
q1a11/max_fs = 775
q1a11/max_ss = 1109
-
q1a12/min_fs = 388
q1a12/min_ss = 1110
q1a12/max_fs = 581
q1a12/max_ss = 1294
-
q1a13/min_fs = 582
q1a13/min_ss = 1110
q1a13/max_fs = 775
q1a13/max_ss = 1294
-
q1a14/min_fs = 388
q1a14/min_ss = 1295
q1a14/max_fs = 581
q1a14/max_ss = 1479
-
q1a15/min_fs = 582
q1a15/min_ss = 1295
q1a15/max_fs = 775
q1a15/max_ss = 1479
-
q2a0/min_fs = 776
q2a0/min_ss = 0
q2a0/max_fs = 969
q2a0/max_ss = 184
-
q2a1/min_fs = 970
q2a1/min_ss = 0
q2a1/max_fs = 1163
q2a1/max_ss = 184
-
q2a2/min_fs = 776
q2a2/min_ss = 185
q2a2/max_fs = 969
q2a2/max_ss = 369
-
q2a3/min_fs = 970
q2a3/min_ss = 185
q2a3/max_fs = 1163
q2a3/max_ss = 369
-
q2a4/min_fs = 776
q2a4/min_ss = 370
q2a4/max_fs = 969
q2a4/max_ss = 554
-
q2a5/min_fs = 970
q2a5/min_ss = 370
q2a5/max_fs = 1163
q2a5/max_ss = 554
-
q2a6/min_fs = 776
q2a6/min_ss = 555
q2a6/max_fs = 969
q2a6/max_ss = 739
-
q2a7/min_fs = 970
q2a7/min_ss = 555
q2a7/max_fs = 1163
q2a7/max_ss = 739
-
q2a8/min_fs = 776
q2a8/min_ss = 740
q2a8/max_fs = 969
q2a8/max_ss = 924
-
q2a9/min_fs = 970
q2a9/min_ss = 740
q2a9/max_fs = 1163
q2a9/max_ss = 924
-
q2a10/min_fs = 776
q2a10/min_ss = 925
q2a10/max_fs = 969
q2a10/max_ss = 1109
-
q2a11/min_fs = 970
q2a11/min_ss = 925
q2a11/max_fs = 1163
q2a11/max_ss = 1109
-
q2a12/min_fs = 776
q2a12/min_ss = 1110
q2a12/max_fs = 969
q2a12/max_ss = 1294
-
q2a13/min_fs = 970
q2a13/min_ss = 1110
q2a13/max_fs = 1163
q2a13/max_ss = 1294
-
q2a14/min_fs = 776
q2a14/min_ss = 1295
q2a14/max_fs = 969
q2a14/max_ss = 1479
-
q2a15/min_fs = 970
q2a15/min_ss = 1295
q2a15/max_fs = 1163
q2a15/max_ss = 1479
-
q3a0/min_fs = 1164
q3a0/min_ss = 0
q3a0/max_fs = 1357
q3a0/max_ss = 184
-
q3a1/min_fs = 1358
q3a1/min_ss = 0
q3a1/max_fs = 1551
q3a1/max_ss = 184
-
q3a2/min_fs = 1164
q3a2/min_ss = 185
q3a2/max_fs = 1357
q3a2/max_ss = 369
-
q3a3/min_fs = 1358
q3a3/min_ss = 185
q3a3/max_fs = 1551
q3a3/max_ss = 369
-
q3a4/min_fs = 1164
q3a4/min_ss = 370
q3a4/max_fs = 1357
q3a4/max_ss = 554
-
q3a5/min_fs = 1358
q3a5/min_ss = 370
q3a5/max_fs = 1551
q3a5/max_ss = 554
-
q3a6/min_fs = 1164
q3a6/min_ss = 555
q3a6/max_fs = 1357
q3a6/max_ss = 739
-
q3a7/min_fs = 1358
q3a7/min_ss = 555
q3a7/max_fs = 1551
q3a7/max_ss = 739
-
q3a8/min_fs = 1164
q3a8/min_ss = 740
q3a8/max_fs = 1357
q3a8/max_ss = 924
-
q3a9/min_fs = 1358
q3a9/min_ss = 740
q3a9/max_fs = 1551
q3a9/max_ss = 924
-
q3a10/min_fs = 1164
q3a10/min_ss = 925
q3a10/max_fs = 1357
q3a10/max_ss = 1109
-
q3a11/min_fs = 1358
q3a11/min_ss = 925
q3a11/max_fs = 1551
q3a11/max_ss = 1109
-
q3a12/min_fs = 1164
q3a12/min_ss = 1110
q3a12/max_fs = 1357
q3a12/max_ss = 1294
-
q3a13/min_fs = 1358
q3a13/min_ss = 1110
q3a13/max_fs = 1551
q3a13/max_ss = 1294
-
q3a14/min_fs = 1164
q3a14/min_ss = 1295
q3a14/max_fs = 1357
q3a14/max_ss = 1479
-
q3a15/min_fs = 1358
q3a15/min_ss = 1295
q3a15/max_fs = 1551
q3a15/max_ss = 1479"""


