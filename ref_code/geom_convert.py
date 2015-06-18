#!/usr/local/bin/env python

'''
This python code offers a handy way to generate CrystFEL and Cheetah formatted geometry files using Mikahil's 
metrology..
@author Shibom Basu..
'''

import sys, os
import numpy as np
import argparse
import math

#First look for options to find out what user wants..
#------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("-e", "--exp", type=str, 
                        help='Experiment number, e.g. cxii6315')
parser.add_argument("-f", "--form", type=str,
                        help="CrystFEL or Cheetah")
parser.add_argument("-o", "--output", type=str, 
                        help="provide filename, for CrystFEL in .geom and for Cheetah in .h5")

args = parser.parse_args()
#------------------------------------------------------------------

#set some necessary constants and parameters to make geometry file..
#------------------------------------------------------------------      
pix_size = 110 # um unit
pix_in_meter = 110e-6 

area_fs = 194; area_ss = 185;
Lside = 388; Sside = 185;

n_asics = 16; # each quadrant has 16 asics or 8 cspad2x1 (the way Mikhail reads it)
n_quads = 4; # cspad is made of 4 quadrants.

hflip = -1

offset_angle = 0
addoffset = 0  #this is a boolean type value. use if needed to apply quadrant offsets..

quadXoffset = [-15,-15,-15,-20];
quadYoffset = [20,20,33,20];
#------------------------------------------------------------------
'''
Now read data from the Mikhail's meterology provided under experimental space
'''
#------------------------------------------------------------------------
if args.exp is not None:
   metrology = "/reg/d/psdm/cxi/"+args.exp+"/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0/geometry/0-end.data"
   print "0-end.data file checked.. \n"
else:
   sys.exit("please provide correct path for 0-end.data file..\n")


cspad = np.loadtxt(metrology, comments='#', usecols = (4,5,6,7,10,11,12))

cspad = cspad[0:32,:]
'''
CrystFEL works in pixel but Mikhail provides geom stuffs in micrometer. So we need to convert..
'''
asics_X = cspad[:,0]/pix_size
asics_Y = cspad[:,1]/pix_size
asics_Z = cspad[:,2]/pix_size

asics_rot = cspad[:,3]; asics_tiltZ = cspad[:,4];
#----------------------------------------------------------------------------------
'''
CrystFEL works with 64 asics total. CrystFEL defines cspad detector as  4 "quadrant" objects,
each of  which contains 16 asics. So, we will generate 4 quadrant objects and calculate required stuffs in 
fast-scan and slow-scan coordinates. Finally return a list of 4 objects.. 
'''


class quadrants(object):
      def __init__(self, nquads, nareas, quadXoffset, quadYoffset, rotation):
            self.nquads = nquads
            self.nareas = nareas
            self.quadXoffset = quadXoffset
            self.quadYoffset = quadYoffset
            self.rotation = rotation
            self.min_fs = []; self.max_fs = [];
            self.min_ss = []; self.max_ss = [];
            self.fs = []; self.ss = [];
            self.cornerX = []; self.cornerY = [];


def rotate2D(v, alpha):
    u = v;
    v[0] = u[0]*math.cos(alpha) - u[1]*math.sin(alpha)
    v[1] = u[0]*math.sin(alpha) + u[1]*math.cos(alpha)
    return v

def pair_num_in_asic(asic_num):
    return round(asic_num/2)

def index_in_pair_asic(asic_num):
    if asic_num % 2 == 0:
       return 0
    else:
       return 1

def crystfel_geom(filename):
       
       
       panels = []; # a list which will contain 4 instances of "quadrant" objects
       ofh = open(filename, 'w')

       ofh.write("photon_energy = /LCLS/photon_energy_ev \n")
       ofh.write("clen = /LCLS/detector0-EncoderValue \n")
       ofh.write("adu_per_ev = 0.00338 \n")
       ofh.write("coffset = 0.587 \n")
       ofh.write("\n")
       ofh.write("; above mentioned fields are often modified for the best results by the crystfel users \n")
       ofh.write("\n\n")


       for q in range(0,n_quads):
          quad = quadrants(q, n_asics, quadXoffset[q], quadYoffset[q], (offset_angle+90 - 90*q))
          panels.append(quad)
    
          for i in range(0,n_asics):
              min_fs = q*2*area_fs + index_in_pair_asic(i)*area_fs
              panels[q].min_fs.append(min_fs)
              min_ss = area_ss*pair_num_in_asic(i)
              panels[q].min_ss.append(min_ss)

              ofh.write("q%ia%i/min_fs = %i\n" %(q,i,min_fs))
              ofh.write("q%ia%i/min_ss = %i\n" %(q,i,min_ss))

              max_fs = q*2*area_fs + index_in_pair_asic(i)*area_fs + area_fs -1
              panels[q].max_fs.append(max_fs)
              max_ss = area_ss*pair_num_in_asic(i) + area_ss -1
              panels[q].max_ss.append(max_ss)

              ofh.write("q%ia%i/max_fs = %i\n" %(q,i,max_fs))
              ofh.write("q%ia%i/max_fs = %i\n" %(q,i,max_ss))

              fs_direction = [0,1]; ss_direction = [1,0];
              angle = asics_rot[pair_num_in_asic(i)] + asics_tiltZ[pair_num_in_asic(i)];
              fs_direct = rotate2D(fs_direction, math.radians(angle));
              ss_direct = rotate2D(ss_direction, math.radians(angle));

              fs_direct = rotate2D(fs_direct, math.radians(quad.rotation))
              ss_direct = rotate2D(ss_direct, math.radians(quad.rotation))

              # flip vertically
              for j in range(0,2): 
                 fs_direct[j] = hflip*fs_direct[j]
                 ss_direct[j] = hflip*ss_direct[j]
        
              panels[q].fs.append(tuple(fs_direct))
              panels[q].ss.append(tuple(ss_direct))

              ofh.write("q%ia%i/fs = %fx + %fy \n" %(q,i,fs_direct[0],fs_direct[1]))
              ofh.write("q%ia%i/ss = %fx + %fy \n" %(q,i,ss_direct[0],ss_direct[1]))
              
              ofh.write("q%ia%i/res = %f \n" %(q,i,9090.91)) #res is a conversion factor that crystfel liks to have

              #first corner of each asic in each quadrant

              center_quad = [0,1]
              center_quad[0] = asics_X[pair_num_in_asic(i)]
              center_quad[1] = asics_Y[pair_num_in_asic(i)]

              corner1_pair = []
              corner1_pair.append(-Sside/2); corner1_pair.append(-Lside/2);
              corner2_pair = []
              corner2_pair.append(-Sside/2); corner2_pair.append(2);

              corner1_pair = rotate2D(corner1_pair, math.radians(angle))
              corner2_pair = rotate2D(corner2_pair, math.radians(angle))

              corner_asic = []
              # corner_asic contains corner_X and corner_Y both for each asic. 
              # first insert corner_X values into corner_asic list

              temp = quad.quadXoffset*addoffset + center_quad[0]
              if index_in_pair_asic(i) == 0:
                 temp += corner1_pair[0]
              else:
                 temp += corner2_pair[0]
              corner_asic.append(temp)

              # now insert corner_Y value in corner_asic list

              temp = quad.quadYoffset*addoffset + center_quad[1]
              if index_in_pair_asic(i) == 0:
                 temp += corner1_pair[1]
              else:
                 temp += corner2_pair[1]
              corner_asic.append(temp)
        
              corner_asic = rotate2D(corner_asic, math.radians(quad.rotation))

              #corner_X
              panels[q].cornerX.append(hflip*corner_asic[0])
              ofh.write("q%ia%i/corner_x = %f \n" %(q,i,(hflip*corner_asic[0])))
              #corner_Y
              panels[q].cornerY.append(corner_asic[1])
              ofh.write("q%ia%i/corner_y = %f \n" %(q,i,corner_asic[1]))
              ofh.write("q%ia%i/badrow_direction = - \n" %(q,i))
              ofh.write("q%ia%i/no_index = 0 \n" %(q,i))

              ofh.write("\n\n")
       ofh.close()
       return panels

def cheetah_geom(filename):
    geoms = GeometryAccess(meterology, 0377)
    '''
    cheetah and crystfel data shape is (1480,1552), i.e., reshape(4,8,185,388)
    '''
    data_shape = (1480,1552)
    X, Y = geoms.get_pixel_coord_indexes()
    X = X.reshape(data_shape)
    Y = Y.reshape(data_shape)
    '''
    Cheetah likes to consider coordinates in physical units i.e., meters instead of pixels..
    '''
    X = np.multiply(X,110e-6);
    Y = np.multiply(Y,110e-6);

    '''
    Cheetah doesn't use Z direction. so, just have it for the sake of having
    Cheetah doesn't need fancy geom file because its just a hitfinding algorithm, 
    where it looks for peak intensities rather than peak location
    '''

    Z = np.zeros(data_shape)
    '''
    Write out geometry file in cheetah format..
    '''
    ofh = h5py.File(filename, 'w')
    dst1 = ofh.create_dataset("x", data=X)#dtype=float32)
    
    dst2 = ofh.create_dataset("y", data=Y)# dtype=float32)
    
    dst3 = ofh.create_dataset("z", data=Z)# dtype=float32)
    
    ofh.close() 

#Lets call the functions and write out geometry files..

if args.form is None:
   sys.exit("tell me the format. \n")

if args.form == 'CrystFEL':
   if args.output is not None:
      filename = args.output
   else:
      filename = "crystfel.geom" 
  
   panel =  crystfel_geom(filename) 

if args.form == 'Cheetah':
   import h5py
   from PSCalib.GeometryAccess import *
   if args.output is not None:
      filename = args.output
   else:
      filename = 'cheetah.h5'
   
   cheetah_geom(filename)


       





