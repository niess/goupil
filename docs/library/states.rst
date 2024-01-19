.. _states:

`states`_
=========

This utility function instantiates a structured
:external:py:class:`numpy.ndarray` of Monte Carlo states that represent photons.


Syntax
------

.. py:function:: states(shape, **kwargs)

   The *shape* argument determines the shape of the
   :external:py:class:`numpy.ndarray`. Optional arguments can be provided to
   initialize specific fields. The available fields are listed below.

   .. list-table::

      * - :python:`"direction"`
        - Current direction (unit vector), using Cartesian coordinates.

      * - :python:`"energy"`
        - Current energy, in MeV.

      * - :python:`"length"`
        - Trajectory path length, in cm.

      * - :python:`"position"`
        - Current position, in cm, using Cartesian coordinates.

      * - :python:`"weight"`
        - Trajectory Monte Carlo weight.


Examples
--------

.. _examples:

Instantiate a size 3 array of photons with default values.

>>> goupil.states(3)
array([...], dtype=...)

As before, but with a specified photon energy of :python:`0.5` MeV for all
states.

>>> goupil.states(3, energy=0.5)
array([...], dtype=...)
