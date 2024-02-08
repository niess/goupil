.. _TopographyMap:

`TopographyMap`_
================

This class wraps a regular grid of elevation data, :math:`z_{ij}`, as function
of :math:`(x_j, y_i)` coordinates at grid nodes. Elevation data are abstracted
as a continuous function

.. math::

   z = f(x, y)

over the grid support by using a bilinear interpolation.

.. warning::

   :doc:`topography_map` objects use indexing order. That is, coordinate
   :python:`(x[j], y[i])` corresponds to elevation :python:`z[i,j]`. This
   differs from reading order used when storing topography data as images, where
   the :math:`y`-axis is inverted.


Constructor
-----------

.. py:class:: TopographyMap(x, y, z=None, shape=None)

   The `x` and `y` arguments specify the support of the grid along :math:`x` and
   :math:`y` coordinates. Without additional arguments, an elevation of
   :python:`0` is considered. For instance,

   >>> zero = goupil.TopographyMap([-1, 1], [-1, 1])

   represents zero elevation (:math:`f = 0`) over the support
   :math:`[-1,1]\times[-1,1]` cm\ :sup:`2`. Optionaly, the `shape` of the
   elevation grid can be specified as a length-2 sequence :math:`(n_y, n_x)`,
   e.g. as

   >>> dem = goupil.TopographyMap([-1, 1], [-2, 2], shape=(21, 41))

   creates an elevation grid with 41 nodes along the :math:`x`-axis and 21 along
   the :math:`y`-axis. A `z` argument can also be provided, as a
   :external:py:class:`float` or :external:py:class:`numpy.ndarray`, in order to
   initialise the elevation grid.


Attributes
----------

.. py:attribute:: TopographyMap.x
   :type: numpy.ndarray

   The grid of :math:`x_j` values.

.. py:attribute:: TopographyMap.y
   :type: numpy.ndarray

   The grid of :math:`y_i` values.

.. py:attribute:: TopographyMap.z
   :type: numpy.ndarray

   The grid of :math:`z_{ij}` values. These values can be modified in-place,
   e.g. as

   >>> dem.z[0,0] = 1.0 # cm

   However, the grid cannot be reshaped.

.. py:attribute:: TopographyMap.box
   :type: tuple

   The map support (bounding-box) along :math:`x` and :math:`y`-axis.


Methods
-------

.. py:method:: TopographyMap.__call__(x, y, grid=None)

   Returns interpolated elevation values at :math:`(x, y)` coordinates. The `x`
   and `y` arguments can be :external:py:class:`float` or
   :external:py:class:`numpy.ndarray` with consistent sizes. If `grid` is
   :python:`True`, then a 2D grid of elevation values is returned over the outer
   product of `x` and `y`. For instance,

   .. doctest::
      :hide:

      >>> xmin, xmax, ymin, ymax = -1, 1, -1, 1

   >>> z = dem(
   ...     numpy.linspace(xmin, xmax, 101),
   ...     numpy.linspace(ymin, ymax, 201),
   ...     grid=True
   ... )

   returns a :math:`201 \times 101` :external:py:class:`numpy.ndarray` of
   elevation values computed over the grid delimited by :math:`[x_\text{min},
   x_\text{max}]\times[y_\text{min}, y_\text{max}]`.
