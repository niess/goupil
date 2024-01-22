.. _RandomStream:

`RandomStream`_
===============

This class represents a pseudo-random stream, such as the one consumed by the
Monte Carlo transport. The stream can be shifted forward or backward using the
:py:attr:`index <RandomStream.index>` attribute. By default, the stream is
seeded using the system entropy.


Constructor
-----------

.. py:class:: RandomStream(seed=None)

   The *seed* argument should be a 128-bit integer. If no *seed* is provided,
   a random value will be chosen using the system entropy.


Attributes
----------

.. note::

   The attributes below are both readable and writable. For instance, decreasing
   the *index* will rewind the pseudo-random stream.

.. py:attribute:: RandomStream.index
   :type: int

   The current index of the pseudo random stream.

.. py:attribute:: RandomStream.seed
   :type: int

   The initial seed of the pseudo random stream.

   .. note::

      Modifying the seed value resets the stream position (index) to
      :python:`0`. If a :python:`None` seed value is set, a random value will be
      selected from the system entropy.


Methods
-------

.. py:method:: RandomStream.normal(shape=None)

   Returns pseudo-random numbers that are distributed according to the Normal
   distribution. If the *shape* argument is not specified, the function will
   return a single :external:py:class:`float`. Otherwise, it will return a
   :external:py:class:`numpy.ndarray` of pseudo-random values with the specified
   shape.

.. py:method:: RandomStream.uniform01(shape=None)

   Returns pseudo-random numbers that are uniformly distributed over
   :math:`(0,1)`. As above, depending on the *shape* argument, the function will
   return a single :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray`.

