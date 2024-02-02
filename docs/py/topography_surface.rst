.. _TopographySurface:

`TopographySurface`_
====================

This class represents a topographic surface as a collection of nested
:doc:`topography_map` objects. The resolution is expected to decrease as one
moves from the most inner map (index :python:`0`) to the most outer one,
providing a Level of Detail (LoD).


Constructor
-----------

.. py:class:: TopographySurface(*maps, offset=None)

   The `maps` arguments must be :doc:`topography_map` objects. The first map
   provided is the innermost. Additionally, a global offset can be specified for
   all elevation values.


Attributes
----------

.. py:attribute:: TopographySurface.maps
   :type: Tuple[TopographyMap]

   The sequence of :doc:`topography_map` objects describing this surface.

.. py:attribute:: TopographySurface.offset
   :type: float

   Global offset applied to all elevation values.

Methods
-------

.. py:method:: TopographySurface.__call__(x, y, grid=None)

   Returns the elevation values at coordinates :math:`(x, y)`, including the
   surface global `offset`. The `x` and `y` arguments can be
   :external:py:class:`float` or :external:py:class:`numpy.ndarray` with
   consistent sizes. If `grid` is set to :python:`True`, elevation values are
   computed over a grid that corresponds to the outer product of `x` and `y`,
   similar to the :py:meth:`TopographyMap.__call__` method.
