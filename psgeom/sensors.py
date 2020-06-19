"""
sensors.py

--- THE PSANA SENSOR CONVENTION

In the "psana" convention, an unrotated reference frame has the following
notations and conventions:

  READOUT   SPACE   AXIS SLICING  |
  -------   -----   ------------  | for `xyz` we assume the slicing is
  slow       -y      xyz[*,:,1]   | (slow, fast, x/y/z/)
  fast       +x      xyz[:,*,0]   |

  so, for example, a JUNGFRAU segment looks like this:

         fast -->
     [0,0]
        -------------------------------------------------
     s  |           |           |           |           |
     l  | 256 x 256 | 256 x 256 | 256 x 256 | 256 x 256 | 
     o  |           |           |           |           | 
     w  -------------------------------------------------  
     |  |           |           |           |           | 
        | 256 x 256 | 256 x 256 | 256 x 256 | 256 x 256 | 
        |           |           |           |           | 
        -------------------------------------------------  
        
   +y ^    
      |
      ---> + x      (+ z out of plane)

"""

import abc
import numpy as np

from psgeom import moveable


# ---- abstract sensor class  --------------------------------------------------


class classproperty(object):
    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)


class SensorElement(moveable.MoveableObject):
    """
    Abstract base class specifying a SensorElement. These elements are the
    actual detecting units of the camera (e.g. a two-by-one for a CSPAD). Many
    such elements usually form a camera.
    
    These objects specify the location of pixels within a sensing element and
    a few other basic facts about them.
    """

    @abc.abstractproperty
    def num_pixels(self):
        return

    @property
    def pixel_shape(self):
        return self._pixel_shape

    @abc.abstractproperty
    def untransformed_xyz(self):
        pass

    @property
    def xyz(self):
        uxyz = self.untransformed_xyz
        T = self.global_transform
        return self._evaluate_transform(T, uxyz)

    @property
    def id_num(self):
        return self._id

    @classmethod
    def from_type(
        cls, type_name, id_num=None, parent=None, rotation_angles=None, translation=None
    ):
        raise NotImplementedError("from_type method not implemented")
        return


# ---- generic sensor implementations  ----------------------------------------


class Gap:
    """
    A simple helper class to define a gap in a SensorArray
    """

    def __init__(self, size, location, axis):
        """
        Parameters
        ----------
        size : float
            The gap size, in pixel units
        location : int
            The location of the gap in the array (integer number of pixels)
        axis : str
            Either 'slow' or 'fast'.
        """
        self.size = float(size)
        self.location = location
        if axis in ["slow", "fast"]:
            self.axis = axis
        else:
            raise ValueError("`axis` must be either `slow` or `fast`")
        return

    @property
    def slc(self):
        """
        This is the slice object that dictates where in the xyz array to
        place the gap
        
        Example
        -------
        xyz[gap.slc] += gap.size
        """
        if self.axis == "fast":
            slc = np.s_[:, self.location :, 0]  # fast/x
        elif self.axis == "slow":
            slc = np.s_[self.location :, :, 1]  # slow/y
        return slc

    @property
    def signed_size(self):
        """
        ... TODO ...
        """
        sign = -1.0 if self.axis == "slow" else 1.0
        return sign * self.size


class PixelArraySensor(SensorElement):
    """
    The PixelArraySensor is an implementation of a rectangular array sensor.
    Likely most cameras can be generated from instances of this element.
    """

    def __init__(
        self,
        shape,
        pixel_shape,
        type_name="None",
        id_num=0,
        parent=None,
        rotation_angles=np.array([0.0, 0.0, 0.0]),
        translation=np.array([0.0, 0.0, 0.0]),
    ):
        """
        Create a PixelArraySensor.
        
        Parameters
        ----------
        shape : tuple of ints
            The shape of the rectangular array of pixels. The ordering of the
            shape is (slow, fast), where fast indicates the direction that is
            most rapidly scanned across when mapping intensities stored in
            a linear memory array onto the camera.
            
        pixel_shape : tuple of floats
            The size of the rectangular pixels that make up the detector, also
            in (slow, fast) directions. Note that units in psgeom are arbitrary,
            but that you need to be consistent!
        
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "QUAD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : PixelArraySensor
            The sensor element.
        """

        self._type_name = type_name
        self._id = id_num

        self.set_parent(parent)

        self._rotation_angles = rotation_angles
        self._translation = translation

        self.shape = tuple(shape)
        self._pixel_shape = np.array(pixel_shape)

        self.gaps = []

        return

    def add_gap(self, size, location, axis):
        """
        Add a gap to the sensor definition.
        
        Parameters
        ----------
        size : float
            The gap size, in pixel units
        location : int
            The location of the gap in the array (integer number of pixels)
        axis : str
            Either 'slow' or 'fast'.
        """
        gp = Gap(size, location, axis)
        self.gaps.append(gp)
        return

    @property
    def num_pixels(self):
        return np.product(self.shape)

    @property
    def dimensions(self):
        """
        Returns
        -------
        dims : 2-tuple of floats
            The lengths along the (slow, fast) scan directions in pixel units.
        """
        slow_gaps_size = sum([g.size for g in self._slow_gaps])
        fast_gaps_size = sum([g.size for g in self._fast_gaps])
        return (
            (self.shape[0] + slow_gaps_size) * self._pixel_shape[0],
            (self.shape[1] + fast_gaps_size) * self._pixel_shape[1],
        )

    @property
    def num_gaps(self):
        return len(self.gaps)

    @property
    def _slow_gaps(self):
        """
        Returns a list of gaps that split along the slow axis, in rev order of
        where they occur along the sensor.
        """
        sgs = [g for g in self.gaps if g.axis == "slow"]
        sgs.sort(key=lambda g: g.location, reverse=True)
        return sgs

    @property
    def _fast_gaps(self):
        """
        Returns a list of gaps that split along the fast axis, in rev order of
        where they occur along the sensor.
        """
        fgs = [g for g in self.gaps if g.axis == "fast"]
        fgs.sort(key=lambda g: g.location, reverse=True)
        return fgs

    def trans_bg_to_sensor(self, bg_data):
        """
        Convert a data array shaped for a basisgrid to one for a sensor,
        which wants to combine segments that have gaps between them. For example,
        for a JUNGFRAU 1M segment, this function will reshape data.
        
          (8,256,256) --> (n,512,1024) 
        
        Parameters
        ----------
        bg_data : np.ndarray
            The data in "basisgrid" format
        
        Returns
        -------
        sensor_data : np.ndarray
            The data in "sensor" format
        
        See Also
        --------
        trans_sensor_to_bg()
            The reverse operation.
        """

        if hasattr(bg_data, "shape"):
            if len(bg_data.shape) == 2:
                bg_data = bg_data.reshape(1, *bg_data.shape)

        s_splits = [sg.location for sg in self._slow_gaps[::-1]]
        f_splits = [fg.location for fg in self._fast_gaps[::-1]]

        n_s = len(s_splits)
        n_f = len(f_splits)

        panels_per_sensor = max(1, 2 * n_s * n_f)

        if len(bg_data) % panels_per_sensor != 0:
            raise ValueError(
                "`bg_data` has %d panels expected %d"
                "" % (len(bg_data), panels_per_sensor)
            )

        ss_data = np.split(np.array(bg_data), n_s + 1, axis=0)
        f_combined = [np.concatenate(e, axis=1) for e in ss_data]
        s_combined = np.concatenate(f_combined, axis=0)
        sensor_data = np.array(s_combined)

        assert sensor_data.shape[-2:] == self.shape, (
            sensor_data.shape[-2:],
            self.shape,
        )

        return np.squeeze(sensor_data)

    def trans_sensor_to_bg(self, sensor_data):
        """
        Convert a data array shaped for a sensor to one for a basisgrid,
        which wants to split segments that have gaps between them. For example,
        for a JUNGFRAU 1M segment, this function will reshape data.
        
          (n,512,1024) --> (8n,256,256)
        
        Parameters
        ----------
        sensor_data : np.ndarray
            The data in "basisgrid" format
        
        Returns
        -------
        bg_data : np.ndarray
            The data in "sensor" format
        
        See Also
        --------
        trans_bg_to_sensor()
            The reverse operation.
        """

        s_splits = [sg.location for sg in self._slow_gaps[::-1]]
        f_splits = [fg.location for fg in self._fast_gaps[::-1]]

        if not sensor_data.shape[-2:] == self.shape:
            raise ValueError(
                "passed data is wrong shape, got %s, expected "
                "%s" % (str(sensor_data.shape, self.shape))
            )

        if len(sensor_data.shape) == 2:
            n_copies = 1
            sensor_data = sensor_data.reshape(1, *sensor_data.shape)
        elif len(sensor_data.shape) == 3:
            n_copies = sensor_data.shape[0]
        else:
            raise ValueError("`sensor_data` must be 2D or 3D array")

        bg_data = []
        for c in range(n_copies):
            s_split_data = np.split(sensor_data[c], s_splits, axis=0)
            for ss in s_split_data:
                bg_data.extend(np.split(ss, f_splits, axis=1))

        shapes = [e.shape for e in bg_data]
        if len(set(shapes)) == 1:
            bg_data = np.array(bg_data)

        return bg_data

    @property
    def untransformed_xyz(self):
        """
        Return the xyz coordinates of the element in the reference frame, that
        is before any translation/rotation operations have been applied.
        """

        xy = np.mgrid[0.0 : float(self.shape[1]), 0.0 : float(self.shape[0])].T
        xy[:, :, :] = xy[::-1, :, :]  # psana convention

        xy[:, :, 0] *= self.pixel_shape[0]
        xy[:, :, 1] *= self.pixel_shape[1]

        # add the z dimension (just flat)
        z = np.zeros([self.shape[0], self.shape[1], 1])
        xyz = np.concatenate([xy, z], axis=-1)

        # add any gaps
        for gap in self.gaps:
            xyz[gap.slc] += gap.signed_size * self.pixel_shape[gap.slc[2]]

        # and, finally, for some reason M [psana] measures rotations from the
        # center of the 2x1 but the corner of the quad. So we center the
        # sensor elements
        xyz[:, :, 0] -= np.mean(xyz[:, :, 0])
        xyz[:, :, 1] -= np.mean(xyz[:, :, 1])

        return xyz

    @property
    def psf(self):
        """
        Return basis grids for this object.
        
        Returns
        -------
        p : np.ndarray
            A 3-vector pointing from the interaction site to the first pixel
            read out from memory for this element.
            
        s : np.ndarray
            A 3-vector pointing along the slow scan direction. The size of the
            vector is the size of the pixel in this direction.
            
        f : np.ndarray
            A 3-vector pointing along the slow scan direction. The size of the
            vector is the size of the pixel in this direction.
        
        shp : tuple
            A 2-tuple of the sensor shape
        """

        xyz = self.xyz

        p = xyz[0, 0, :]
        s = xyz[1, 0, :] - p
        f = xyz[0, 1, :] - p

        grids = []
        slow_split_grids = []

        # >>> slow scan
        # for each gap along the slow axis create a new grid

        curr_shp = self.shape  # track how much is left to divide up
        for ig, gap in enumerate(self._slow_gaps):  # gaps come in rev order

            # we need to count the spacing for all the "upstream" gaps
            tot_gap_size = sum([gx.size for gx in self._slow_gaps[ig:]])

            new_p = p + s * (gap.location + tot_gap_size)
            new_shp = (curr_shp[0] - gap.location, curr_shp[1])
            curr_shp = (gap.location, curr_shp[1])

            slow_split_grids.append([new_p, s, f, new_shp])

        # add the remaining (first) panel
        slow_split_grids.append([p, s, f, curr_shp])

        # >>> fast scan
        # then, for each grid, split along the fast axis gaps

        for grid in slow_split_grids:

            p = grid[0]  # use shifted value
            curr_shp = grid[3]
            for ig, gap in enumerate(self._fast_gaps):  # gaps come in rev order

                # we need to count the spacing for all the "upstream" gaps
                tot_gap_size = sum([gx.size for gx in self._fast_gaps[ig:]])

                new_p = p + f * (gap.location + tot_gap_size)
                new_shp = (curr_shp[0], curr_shp[1] - gap.location)
                curr_shp = (curr_shp[0], gap.location)

                grids.append([new_p, s, f, new_shp])

            # add the remaining (first) panel
            grids.append([grid[0], grid[1], grid[2], curr_shp])

        assert len(grids) == max(2 * self.num_gaps, 1), (len(grids), self.num_gaps)

        ret = [tuple(g) for g in grids]  # convert to tuples

        # we iterated through the gaps in reverse-position order
        # so reverse the list to get increasing-position ordering
        ret.reverse()
        assert np.all(ret[0][0] == xyz[0, 0, :])  # first p is same as old p

        return ret

    @classmethod
    def from_type(
        cls,
        type_name,
        id_num=0,
        parent=None,
        rotation_angles=np.array([0.0, 0.0, 0.0]),
        translation=np.array([0.0, 0.0, 0.0]),
    ):
        """
        Factory function for automatically identifying
        the sensor based on the `type_name` alone.            
        """
        # see code in translate.py
        return cls(
            type_name=type_name,
            id_num=id_num,
            parent=parent,
            rotation_angles=rotation_angles,
            translation=translation,
        )


class FixedArraySensor(PixelArraySensor):
    """
    Slight modification of the PixelArraySensor that requires a fixed sensor
    shape and pixel size. Should be preferred when these are known.
    """

    def __init__(self, **kwargs):
        super(FixedArraySensor, self).__init__(self.shape, self.pixel_shape, **kwargs)
        return

    @abc.abstractmethod
    def shape(self):
        # implement this as a @classproperty
        return

    @abc.abstractmethod
    def pixel_shape(self):
        # implement this as a @classproperty
        return


# ---- specific sensor implementations  ---------------------------------------


class MtrxV1(PixelArraySensor):
    """
    This class is to ensure backwards compatability for "MTRX:V1" sensor 
    elements.

	These "V1" elements have a different convention as to how the origin
	and slow/fast axes are placed: 
       * the origin is at the first pixel
       * the ss direction is along +x
       * the fs direction is along +y

    The overwritten method "untransformed_xyz" takes care of this.
    """
    
    def __init__(self, shape, pixel_shape, id_num=0, parent=None,
                 rotation_angles=np.array([0.0, 0.0, 0.0]), 
                 translation=np.array([0.0, 0.0, 0.0])):
        """
        Create a Mtrx.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "RAYONIX:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : Mtrx
            The sensor element.
        """
                 
        if shape is None or pixel_shape is None:
            raise RuntimeError('shape or pixel shape not supplied to Mtrx')

        # TJL 4/9/18
        # I am not sure why these lines are necessary
        # but they seem to be to get these attributes set
        # I would have expected the super init method below to take care of it...
        self.shape = shape
        self._pixel_shape = pixel_shape

        super(MtrxV1, self).__init__(shape, pixel_shape, 
                 type_name='shouldbeoverwritten', 
                 id_num=id_num, parent=parent,
                 rotation_angles=rotation_angles, 
                 translation=translation)
                                       
        return


    @property
    def untransformed_xyz(self):
        """
        Return the xyz coordinates of the element in the reference frame, that
        is before any translation/rotation operations have been applied.
        """
                
        # convention that x/row/first-index is the SLOW varying dimension
        # y/column/second-index is FAST and z is perpendicular to the sensor 
        # completing a right handed coordinate system in the untransformed view
        
        xy = np.mgrid[0.0:float(self.shape[0]),0.0:float(self.shape[1])]
        xy = np.rollaxis(xy, 0, start=3)
        
        xy[:,:,0] *= self.pixel_shape[0]
        xy[:,:,1] *= self.pixel_shape[1]
        
        # add the z dimension (just flat)
        z = np.zeros([self.shape[0], self.shape[1], 1])
        xyz = np.concatenate([xy, z], axis=-1)
        
        return xyz   


    @classmethod
    def from_type(cls, type_name,
                  id_num=0, parent=None,
                  rotation_angles=np.array([0.0, 0.0, 0.0]),
                  translation=np.array([0.0, 0.0, 0.0])):
        """
        Factory function for automatically identifying
        the sensor based on the `type_name` alone.            
        """

        if 'MTRX' not in type_name:
            raise ValueError('`type_name` (%s) does not contain "MTRX"'
                             'cannot generate Mtrx object from type_name'
                             ' alone')

        s0, s1, ps0, ps1 = type_name.split(':')[-4:]
        shape = (int(s0), int(s1))
        pixel_shape = (float(ps0), float(ps1))

        return cls(shape, pixel_shape,
                   id_num=id_num, parent=parent,
                   rotation_angles=rotation_angles,
                   translation=translation)

    @property
    def type_name(self):
        return 'MTRX:%d:%d:%d:%d' %(self.shape[0], self.shape[1], 
                                    round(self._pixel_shape[0]),
                                    round(self._pixel_shape[1]))


class Mtrx(PixelArraySensor):
    """
    A specific PixelArraySensor representing a generic rectangular sensor.
    """

    def __init__(
        self,
        shape,
        pixel_shape,
        id_num=0,
        parent=None,
        rotation_angles=np.array([0.0, 0.0, 0.0]),
        translation=np.array([0.0, 0.0, 0.0]),
    ):
        """
        Create a Mtrx.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "RAYONIX:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : Mtrx
            The sensor element.
        """

        if shape is None or pixel_shape is None:
            raise RuntimeError("shape or pixel shape not supplied to Mtrx")

        # TJL 4/9/18
        # I am not sure why these lines are necessary
        # but they seem to be to get these attributes set
        # I would have expected the super init method below to take care of it...
        self.shape = shape
        self._pixel_shape = pixel_shape

        super(Mtrx, self).__init__(
            shape,
            pixel_shape,
            type_name="shouldbeoverwritten",
            id_num=id_num,
            parent=parent,
            rotation_angles=rotation_angles,
            translation=translation,
        )

        return

    @classmethod
    def from_type(
        cls,
        type_name,
        id_num=0,
        parent=None,
        rotation_angles=np.array([0.0, 0.0, 0.0]),
        translation=np.array([0.0, 0.0, 0.0]),
    ):
        """
        Factory function for automatically identifying
        the sensor based on the `type_name` alone.            
        """

        if "MTRX" not in type_name:
            raise ValueError(
                '`type_name` (%s) does not contain "MTRX"'
                "cannot generate Mtrx object from type_name"
                " alone"
            )

        s0, s1, ps0, ps1 = type_name.split(":")[-4:]
        shape = (int(s0), int(s1))
        pixel_shape = (float(ps0), float(ps1))

        return cls(
            shape,
            pixel_shape,
            id_num=id_num,
            parent=parent,
            rotation_angles=rotation_angles,
            translation=translation,
        )

    @property
    def type_name(self):
        return "MTRX:V2:%d:%d:%d:%d" % (
            self.shape[0],
            self.shape[1],
            round(self._pixel_shape[0]),
            round(self._pixel_shape[1]),
        )


class Cspad2x1(FixedArraySensor):
    """
    CSPAD 2x1 panel
    """

    @classproperty
    def shape(self):
        return (185, 388)

    @classproperty
    def pixel_shape(self):
        return np.array([109.92, 109.92])  # microns

    def __init__(self, **kwargs):
        """
        Create a Cspad2x1.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "QUAD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : Cspad2x1
            The sensor element.
        """

        if "type_name" not in kwargs.keys():
            kwargs["type_name"] = "SENS2X1:V1"
        super(Cspad2x1, self).__init__(**kwargs)
        self.add_gap(3.0, 194, "fast")

        return


class PnccdQuad(FixedArraySensor):
    """
    A pnCCD quad.
    """

    @classproperty
    def shape(self):
        return (512, 512)

    @classproperty
    def pixel_shape(self):
        return np.array([75.0, 75.0])  # microns

    def __init__(self, **kwargs):
        """
        Create a PnccdQuad.
        
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "PNCCD:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : PnccdQuad
            The sensor element.
        """

        if "type_name" not in kwargs.keys():
            kwargs["type_name"] = "PNCCD:V1"
        super(PnccdQuad, self).__init__(**kwargs)

        return


class JungfrauSegment(FixedArraySensor):
    """
    A specific PixelArraySensor representing a 1M JUNGFRAU segment.
    """

    @classproperty
    def shape(self):
        return (512, 1024)

    @classproperty
    def pixel_shape(self):
        return np.array([75.0, 75.0])  # microns

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "JUNGFRAU:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.
            
        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : JungfrauSegment
            The sensor element.
        """

        if "type_name" not in kwargs.keys():
            kwargs["type_name"] = "JUNGFRAU:V1"
        super(JungfrauSegment, self).__init__(**kwargs)

        self.add_gap(2.0, 256, "fast")
        self.add_gap(2.0, 512, "fast")
        self.add_gap(2.0, 768, "fast")
        self.add_gap(2.0, 256, "slow")

        return


class Epix10kaSegment(FixedArraySensor):
    """
    A specific PixelArraySensor representing a EPIX10KA segment.

    The EPIX10KA segment has size 352x384 and is composed of 4
    ASICs. Each ASIC is actually 178x192, but the last two row are
    unbound, used for calibration and automatically removed from
    the array provided by psana. The pixels at the edge of the
    ASICs are "big". Pixels between ASICs have a size 100x250 µm.
    The four central pixels of the segment have size 250x250 µm².
    """

    @classproperty
    def shape(self):
        return (352, 384)

    @classproperty
    def pixel_shape(self):
        return np.array([100.0, 100.0])  # microns

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        type_name : str
            Give this detector a descriptive name. Often there might be
            two different instances of CompoundDetector with the same name,
            if they are identical units. E.g., "EPIX10KA:V1".
            
        id_num : int
            The unit should have an index. This is not only a unique identifier
            but helps order elements within the camera tree, which can change
            the way someone wants to map pixel intensities (somewhere else in
            memory) onto the camera geometry.
            
        parent : CompoundDetector
            The parent frame, specified by an instance of CompoundDetector.
            
        rotation_angles : np.ndarray
            Three Cardan angles specifying the local frame rotation operator.
            Argument must be a one-D 3-vector.

        translation : np.ndarray
            The xyz translation of the local frame. Argument must be a one-D 
            3-vector.
            
        Returns
        -------
        self : JungfrauSegment
            The sensor element.
        """

        if "type_name" not in kwargs.keys():
            kwargs["type_name"] = "EPIX10KA:V1"
        super(Epix10kaSegment, self).__init__(**kwargs)

        self.add_gap(3.0, 192, "fast")
        self.add_gap(3.0, 176, "slow")

        return

