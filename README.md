# psgeom
[![Build Status](https://travis-ci.org/slaclab/psgeom.svg?branch=master)](https://travis-ci.org/slaclab/psgeom)

psgeom aims to provide an easy to use code base for common geometrical 
operations during scattering experiments:
* load/save multiple geometry formats
* translate/rotate parts of a detector
* easily compute reciprocal/polar space coordinates
* perform angular integration

User-friendliness is emphasized. The software aims to be general
but not sacrifice simplicity.

Chances are that if you need to compute stuff relating to scattering
geometry for a particular experiment, psgeom has what you need, and
it will be easy to use.

Currently, interfaces exist to geometry formats from:
* psana
* CrystFEL
* cheetah
* DIALS
* LCLS detector group metrologies
* a simple flat text / HDF5 pixel map

TJ Lane <thomas.joseph.lane@gmail.com>

------

## Scripts ##

A lot of users just want to convert geometry files between different formats, or other simple tasks. To assist, psgeom provides a few command line scripts that will be of general interest:

* `geoconv` : convert between geometry file formats
* `geoQ` : quickly get reciprocal space coordinates for a geometry
* `gainmk` and `gainconv` : for CSPAD gain files -- to make a new one and convert it's formatting, respectively

Run any script with a `-h` flag to get more information

------

## Examples ##
A few examples of how to use the code.

### format conversions ###
```python
from psgeom import camera

geom = camera.CompoundAreaCamera.from_psana_file('1-end.data')
geom.to_crystfel_file('my_new.geom')
```

### looking at pixel positions ###
```python
geom = camera.CompoundAreaCamera.from_crystfel_file('my.geom')
print(geom.xyz) # real-space xyz coords
```

### looking at things in basisgrid format ###
By "basisgrid", we mean a the geometry is described as a set of panels; each panel by a vector pointing
to the first pixel to be read from memory, along with two vectors for the slow/fast scan directions.
```python
bg = geom.to_basisgrid()
for g in bg.num_grids:
    print(bg.get_grid(g))
```

### radial averaging ###
```python
from psgeom import bin

xyz = geom.xyz

beam_vector = np.array([0.0, 0.0, 1.0])     # assumed
wavenumber = 2.0 * np.pi / wavelength       # inv A

norm = np.linalg.norm(xyz, axis=-1)
S = xyz / norm[...,None] # unit vector

q_xyz = wavenumber * (S - beam_vector)
q_mag = np.linalg.norm(q_xyz, axis=-1)

radavg = bin.Averager(q_mag, mask, n_bins=500)

# data is raw detector data, Iq is radial average in q-coords
Iq1 = radavg(data1) 
Iq2 = radavg(data2) 
...
```

### an easier way to compute reciprocal coords ###
```python
from psgeom import camera
from psgeom import reciprocal

geom = camera.CompoundAreaCamera.from_crystfel_file('my.geom')
d = reciprocal.Geometry(geom)

d.xyz        # real space cart
d.polar      # real space polar
d.reciprocal # reciprocal cart
d.recpolar   # reciprocal polar coords

# compute interpolated values at specific intersection points
pix, intersect = d.compute_intersections(interp_vectors, 0) # 0 --> grid_index
```


-------

## How it works ##

In the scattering community, there are two principle paradigms in use for representing complex detector/experimental geometries: the "basisgrid" paradigm, and the "elemental" paradigm. `psgeom` implements *both*, which means it is both extremely powerful and capable of converting between formats that use either paradigm.

The "basisgrid" paradigm represents a geometry as a set of 2d planar sensors in 3d space. These sensors are represented by a position (`p`) vector that points to the first pixel of that panel to be read from memory. Two additional vectors correspond to the slow (`s`) and fast (`f`) scan directions, showing how to map a 2d array of data onto the sensor geometry. Used by, for example, crystFEL.

The "elemental" paradigm defines sensor components by hand, and then builds a more complicated sensor geometry by translating/rotating these elements with respect to one another, perhaps in a heirarchical fashion. This provides a bit extra power in terms of representing how rigid units move in space, as they are often correlated (e.g. mounted on a hard platform driven by motors in a beamline). Because the elements must be specified in software, this paradigm comes with additional complexity.

In reality all known x-ray and electron detectors to date (as of April 2020) are , composed of planar, 2d arrays of rectangular pixels. `psgeom` reflects this by emphasizing element definitions that correspond to this pattern; such elements are defined by their (1) pixel shape, (2) array size (n-by-m pixels), and (3) any gaps in the sensor surface. This allows for conversion between the "elemental" and "basisgrid" paradigms.

A paper describing these formats and the mathematics behind them is in preparation

-------

### Quick Links ###

Information about the psana geometry:
https://confluence.slac.stanford.edu/display/PSDM/Detector+Geometry

List of high quality geometries generated by users and optical metrologies generated by LCLS:
https://confluence.slac.stanford.edu/display/PSDM/Geometry+History

CrystFEL Geometry:
http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html


