.. _SphereShape:

`SphereShape`_
==============

This class represents a geometric shape in the form of a sphere. It can be used
to define a collection :py:attr:`boundary <TransportSettings.boundary>` for
forwards transported gamma-rays, or as a generation surface in the backward
case.


Constructor
-----------

.. py:class:: SphereShape(radius, center=None)

   This function creates a sphere with specified *radius*, in cm. By default,
   the sphere is centered on the origin. The *centre* argument can be used to
   specify an offset, in cm.

   For example, the following creates a sphere of 1m radius with a 5 cm offset
   along the $z$-axis

   >>> shape = goupil.SphereShape(
   ...     1e2, # cm
   ...     center = (0.0, 0.0, 5.0) # cm
   ... )

Attributes
----------

.. py:attribute:: SphereShape.center
   :type: numpy.ndarray

   The Cartesian coordinates of the sphere center, in cm, as a size-3
   :external:py:class:`numpy.ndarray`.

.. py:attribute:: SphereShape.radius
   :type: numpy.ndarray | None

   The sphere radius, in cm.


Methods
-------

.. py:method:: SphereShape.distance(states, reverse=None) -> numpy.ndarray

   Returns a :external:py:class:`numpy.ndarray` of floats that indicate the
   distance of Monte Carlo *states* to the sphere surface, along their
   respective momentum directions.

   By default, Monte Carlo states are assumed to propagate forwards. To trace
   states backwards instead, set the reverse flag to :python:`False`.

.. py:method:: SphereShape.inside(states) -> numpy.ndarray

   Returns a :external:py:class:`numpy.ndarray` of booleans that indicate
   whether the provided Monte Carlo *states* are inside the sphere or not.

.. py:method:: SphereShape.sample(states, engine=None, side=None, direction=None, weight=None)

   Samples Monte Carlo states over the sphere surface. The *weight* boolean flag
   indicates wether the Monte Carlo states should be weighted by the inverse of
   the sampling PDF or not.

   The *side* arguments refers to the surface side on which positions are
   sampled, within a numerical epsilon. Possible values are are
   :python:`"Inside"` or :python:`"Outside"`.

   The *direction* arguments indicates the orientation of the sampled Monte
   Carlo states with respect to the sphere surface. Possible values are
   :python:`"Ingoing"` or :python:`"Outgoing"`. By default, Monte Carlo states
   are considered as ingoing.

   If a :doc:`transport_engine` is provided (using the *engine* argument), then
   the sampling is configured based one the simulation :py:attr:`mode
   <TransportSettings.mode>`. In a a forward simulation Monte Carlo states are
   generated on the inside by default, while in the backward case, the outside
   is used.

.. note::

   At return from the :py:meth:`sample <SphereShape.sample>` method, the
   positions and directions of Monte Carlo states are modified in-place, along
   with their weights if applicable.

