.. _BoxShape:

`BoxShape`_
===========

This class represents a geometric shape in the form of a box. It can be used to
define a collection :py:attr:`boundary <TransportSettings.boundary>` for
forwards transported gamma-rays, or as a generation surface in the backward
case.


Constructor
-----------

.. py:class:: BoxShape(size, center=None, rotation=None)

   This function creates a box with specified dimensions, in cm. The *size* can
   be given as a :external:py:class:`float` (for a cube) or as a size-3 sequence
   of floats, indicating the width (along *x*, *y*) and height (along *z*) of
   the box. By default, the box is axis-aligned and centered on the origin. The
   *centre* argument can be used to specify an offset, in cm. Additionally, a
   3x3 *rotation* matrix can be specified.

   For example, the following creates a cubic box of 1m width with a 5 cm offset
   along the $z$-axis

   >>> shape = goupil.BoxShape(
   ...     1e2, # cm
   ...     center = (0.0, 0.0, 5.0) # cm
   ... )

Attributes
----------

.. py:attribute:: BoxShape.center
   :type: numpy.ndarray

   The Cartesian coordinates of the box center, in cm, as a size-3
   :external:py:class:`numpy.ndarray`.

.. py:attribute:: BoxShape.rotation
   :type: numpy.ndarray | None

   An optional 3x3 rotation matrix specifying the box orientation, if not
   axis-aligned.

.. py:attribute:: BoxShape.size
   :type: numpy.ndarray

   The box dimensions along the *x*, *y* and *z*-axis (before any rotation), in
   cm.


Methods
-------

.. py:method:: BoxShape.distance(states, /, *, reverse=None) -> numpy.ndarray

   Returns a :external:py:class:`numpy.ndarray` of floats that indicate the
   distance of Monte Carlo *states* to the box surface, along their respective
   momentum directions.

   By default, Monte Carlo states are assumed to propagate forwards. To trace
   states backwards instead, set the reverse flag to :python:`False`.

.. py:method:: BoxShape.inside(states, /) -> numpy.ndarray

   Returns a :external:py:class:`numpy.ndarray` of booleans that indicate
   whether the provided Monte Carlo *states* are inside the box or not.

.. py:method:: BoxShape.sample(states, /, *, engine=None, side=None, direction=None, weight=None)

   Samples Monte Carlo states over the box surface. The *weight* boolean flag
   indicates wether the Monte Carlo states should be weighted by the inverse of
   the sampling PDF or not.

   The *side* arguments refers to the surface side on which positions are
   sampled, within a numerical epsilon. Possible values are
   :python:`"Inside"` or :python:`"Outside"`.

   The *direction* arguments indicates the orientation of the sampled Monte
   Carlo states with respect to the box surface. Possible values are
   :python:`"Ingoing"` or :python:`"Outgoing"`. By default, Monte Carlo states
   are considered as ingoing.

   If a :doc:`transport_engine` is provided (using the *engine* argument), then
   the sampling is configured based one the simulation :py:attr:`mode
   <TransportSettings.mode>`. In a a forward simulation Monte Carlo states are
   generated on the inside by default, while in the backward case, the outside
   is used.

.. note::

   At return from the :py:meth:`sample <BoxShape.sample>` method, the positions
   and directions of Monte Carlo states are modified in-place, along with their
   weights if applicable.

