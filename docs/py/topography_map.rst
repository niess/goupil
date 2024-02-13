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

.. py:class:: TopographyMap(x, y, z=None)

   The `x` and `y` arguments specify the support of the grid along :math:`x` and
   :math:`y` coordinates. Without additional arguments, an elevation of
   :python:`0` is considered. For instance,

   >>> zero = goupil.TopographyMap([-1e5, 1e5], [-1e5, 1e5])

   represents zero elevation (:math:`f = 0`) over the support
   :math:`[-1,1]\times[-1,1]` km\ :sup:`2`. A `z` argument can also be provided,
   as a :external:py:class:`float` or 2D :external:py:class:`numpy.ndarray`, in
   order to initialise the elevation grid.

.. note::

   If a grid node has no data, :external:py:data:`numpy.nan` should be used.


Attributes
----------

.. note::

   :doc:`topography_map` attributes are read-only.

.. py:attribute:: TopographyMap.x
   :type: numpy.ndarray

   The grid of :math:`x_j` values.

.. py:attribute:: TopographyMap.y
   :type: numpy.ndarray

   The grid of :math:`y_i` values.

.. py:attribute:: TopographyMap.z
   :type: numpy.ndarray

   The grid of :math:`z_{ij}` values.

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
      >>> topography = zero

   >>> z = topography(
   ...     numpy.linspace(xmin, xmax, 101),
   ...     numpy.linspace(ymin, ymax, 201),
   ...     grid=True
   ... )

   returns a :math:`201 \times 101` :external:py:class:`numpy.ndarray` of
   elevation values computed over the grid delimited by :math:`[x_\text{min},
   x_\text{max}]\times[y_\text{min}, y_\text{max}]`.
