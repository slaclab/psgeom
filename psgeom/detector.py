#!/usr/bin/env python

"""
Detector Geometries
===================

Briefly, geometries are organized heriarchically: there are `SensorElements`
that compose the leaves of a tree. These represent arbitrary objects that
actually measure e.g. photon intensities of scattered x-rays. An example of
such an element is an ASIC on a Pilatus, or a 2x1 on a CSPAD.

These elements are then "mounted" onto `CompoundDetector`s, which represent
physical units that translate and rotate together. For example, 8 CSPAD 2x1s
are mounted on a single "Quad", that moves as a collective unit. The elements
composing a CompoundDetector` instance can be SensorElements or other 
CompoundDetectors.

A note about how the heirarchical geometry orientation is applied. Each node
in the graph contains a rotation and translation with respect to it's parent.
The translation is applied *first*. So if T is a translation operator, and
R is a rotation,

    x_final = Rx + p

In practice, both are performed using matrices.


Author: TJ Lane <tjlane@slac.stanford.edu>
June 11, 2015



To Fix
------
--> mikhail file z/y/x & euler angle conventions
--> document slow/fast scan x/y conventions
--> think about how to separate sensors out
--> test export methods
--> connect BasisGrid, test read methods

--> fill in legacy.py

--> write "cspad class" to recognize slight variants (4,16,...), flat 
    hierarchies and convert them to cannonical CSPAD form

--> think about the mapping of intensity data onto the detector

To Do
-----
1) implement more general sketch, assemble
2) implement show image on 2d projection
3) think: should mask be included?
5) dynamically expose leaf properties to parents?

"""

import warnings
import numpy as np

from psgeom import moveable
from psgeom import sensors
from psgeom import translate
from psgeom import basisgrid

        
class CompoundDetector(moveable.MoveableObject, moveable.MoveableParent):
    
    def __init__(self, type_name=None, id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):
        
        self._type_name = type_name
        self._id        = id_num
        
        self.set_parent(parent)
        self._children = []
        
        self._rotation_angles = rotation_angles
        self._translation     = translation
        
        return
        
        
    def add_child(self, child):
        """
        Add a child to the compound detector. This can either be a
        `SensorElement` or another `CompoundDetector`.
        
        Parameters
        ----------
        child : SensorElement or CompoundDetector
            The child object to add to this CompundDetector node.
        """
        
        if not (isinstance(child, CompoundDetector) or \
                isinstance(child, sensors.SensorElement)):
                raise TypeError('`child` must be type: SensorElement or '
                                'CompoundDetector')
        
        for c in self.children:
            if c.name == child.name:
                if c is child:
                    raise NameError('Child object already registered with parent!')
                else:
                    raise NameError('Child with name %s already registered with'
                                    ' this parent (%s) -- please change the ID'
                                    ' number to give this object a unique name '
                                    'and re-register it as a child object' % \
                                    (child.name, self.name))
        
        self.children.append(child)
        child._parent = self
        
        return
        
        
    def draw_tree(self):

        print "--- " + str(self.name)
        
        def draw_child_tree(current, depth):
        
            for c in current.children:
                print depth * "    " + "|-- " + str(c.name)
                if hasattr(c, 'children'):
                    draw_child_tree(c, depth + 1)
                    
        draw_child_tree(self, 1)
        
        return
    
        
    @property
    def num_children(self):
        return len(self._children)
    
        
    @property
    def children(self):
        return self._children
        
        
    @property
    def leaves(self):
        
        leaves = []
        
        def add_leaves(node):
            for c in node.children:
                if hasattr(c, 'children'):
                    add_leaves(c)
                else:
                    leaves.append(c)
                    
        add_leaves(self)
                    
        return leaves
        
    
    @property
    def num_pixels(self):
        return np.sum([ c.num_pixels for c in self._children ])
        
        
    @property
    def xyz(self):
        return np.array([ c.xyz for c in self._children ])


    # begin interfaces -----
    
    def to_basisgrid(self):
        
        bg = basisgrid.BasisGrid()
        
        for sensor in self.leaves:
            if not isinstance(sensor, sensors.PixelArraySensor):
                raise TypeError('basisgrid representation is only compatible '
                                'with detectors that are entirely comprised of '
                                'PixelArrayElements')
                               
            p, s, f = sensor.psf 
            bg.add_grid(p, s, f, sensor.shape)
        
        return bg
        
    
    @classmethod
    def from_basisgrid(cls, bg):
        """
        docstring
        """
        
        if not isinstance(bg, basisgrid.BasisGrid):
            raise TypeError('`bg` argument must be instance of BasisGrid,'
                            ' got: %s' % type(bg))
        
        cd = cls(type_name='root_frame', id_num=0, parent=None)
        
        
        for g in range(bg.num_grids):
            
            p, s, f, shape = bg.get_grid(g)
            
            pixel_shape = (np.linalg.norm(s),
                           np.linalg.norm(f))
            
            # to compute the rotation, find the 
            us = s / pixel_shape[0] # unit vector
            uf = f / pixel_shape[1] # unit vector
            n  = np.cross(uf, us)   # tested for orthog. in next fxn
            
            ra = moveable._angles_from_rotated_frame(uf, us, n)

            # translation is just p
            tr = p
            
            pas = sensors.PixelArraySensor(shape, 
                                           pixel_shape, 
                                           type_name='grid_element_%dx%d' % shape, 
                                           id_num=g, 
                                           parent=cd,
                                           rotation_angles=ra, 
                                           translation=tr)
        
        
        return cd
    
    
    def to_psana_file(self, filename, title='geometry'):
        """
        """
        translate.write_psana(self, filename, title)
        return
        
        
    @classmethod
    def from_psana_file(cls, filename):
        """
        """
        return translate.load_psana(cls, filename)
        
    
    def to_crystfel_file(self, filename):
        """
        Write a geometry to disk in CrystFEL format. Note that some fields
        will be written but left blank -- these are fields you probably should
        fill in before performing any computations in CrystFEL, but are 
        information that we have no handle on (e.g. detector gain).

        Thanks to Rick Kirian & Tom White for assistance with this function.

        Parameters
        ----------
        filname : str
            The name of file to write. Will end in '.geom'

        Optional Parameters
        -------------------
        intensity_file_type : str, {'cheetah'}
            The kind of file this geometry file will be used with. Necessary to 
            tell CrystFEL how intensity data map onto the detector
        """
        translate.write_crystfel(self, filename, intensity_file_type='cheetah')
        return
        
        
    @classmethod
    def from_crystfel_file(cls, filename):
        # translate.load_crystfel(cls, filename)
        raise NotImplementedError()
        
        

    

# ---- specific detector implementations ---------------------------------------

class CSPAD(CompoundDetector):
        
        
    def to_cheetah_file(self, filename):
        """
        Convert a CompoundDetector object to a Cheetah h5 pixel map.

        The CompoundDetector must be the a (4,8,185,388) shape representing a
        CSPAD -- cheetah can only understand CSPAD geometries. 

        Parameters
        ----------
        geometry : cspad.CSPad
            The detector geometry to write to disk
    	filename : string
    		The file name for the output pixel map
        """
        translate.write_cheetah(self, filename)
        return


    def from_cheetah_file(filename):
        raise NotImplementedError()




