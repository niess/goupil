.. _TopographySurface:

`TopographySurface`_
====================

This class represents a pseudo-random stream, such as the one consumed by the
Monte Carlo transport. The stream can be shifted forward or backward using the
:py:attr:`index <RandomStream.index>` attribute. By default, the stream is
seeded using the system entropy.


Constructor
-----------

.. py:class:: TopographySurface(*maps, offset=None)

   The *seed* argument should be a 128-bit integer. If no *seed* is provided,
   a random value will be chosen using the system entropy.


Attributes
----------

.. py:attribute:: TopographySurface.offset
   :type: float


Methods
-------

.. py:method:: TopographySurface.__call__(x, y, grid=None)

   Returns pseudo-random numbers that are distributed according to the Normal
   distribution. If the *shape* argument is not specified, the function will
   return a single :external:py:class:`float`. Otherwise, it will return a
   :external:py:class:`numpy.ndarray` of pseudo-random values with the specified
   shape.
