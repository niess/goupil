:class:`DensityGradient`
========================

.. _DensityGradient:

----

This class represents an exponential density gradient, for example, for
modelling the density of the atmosphere. The density varies continuously with
location :math:`\vec{r}` as

.. math::

   \rho(\vec{r}) = \rho_0 \exp\left((\vec{r}-\vec{r}_0)\cdot \frac{\vec{n}}{\lambda}\right),

where :math:`\rho_0`, :math:`\vec{r}_0`, :math:`\vec{n}` and :math:`\lambda` are
configurable model parameters (see the constructor below).


Constructor
-----------

.. py:class:: DensityGradient(density, scale, direction=None, origin=None)

   Creates a density gradient along a given *direction* (:math:`\vec{n}`). The
   *density* argument (:math:`\rho_0`, in :math:`\text{g}/\text{cm}^3`)
   specifies the density value at the *origin* (:math:`\vec{r}_0`, in
   :math:`\text{cm}`). The *scale* argument (:math:`\lambda`, in
   :math:`\text{cm}`) controls the density variation, as detailed above. The
   optional *direction* and *origin* arguments should be 3d-cartesian
   coordinates (e.g. a :external:py:class:`tuple`), refering to the simulation
   frame. If no *direction* is provided, then :python:`(0,0,-1)` is assumed. By
   default, the gradient *origin* is at :python:`(0,0,0)`.

Attributes
----------

.. py:attribute:: DensityGradient.density
   :type: float

   The density value (:math:`\rho_0`) at the gradient :py:attr:`origin
   <DensityGradient.origin>`, in :math:`\text{g}/\text{cm}^3`.

.. py:attribute:: DensityGradient.direction
   :type: numpy.ndarray

   The direction (:math:`\vec{n}`) of the density gradient, in Cartesian
   coordinates.

.. py:attribute:: DensityGradient.origin
   :type: numpy.ndarray

   The origin (:math:`\vec{r}_0`) of the density gradient model, in Cartesian
   coordinates (in :math:`\text{cm}`).

.. py:attribute:: DensityGradient.scale
   :type: float

   The scale parameter (:math:`\lambda`) of the exponential density gradient
   model (in :math:`\text{cm}`).


Methods
-------

.. py:method:: DensityGradient.__call__(position)

   Returns the density value(s) at the requested position(s). The *position*
   argument can be a length 3 sequence or a shape :math:`[\cdots, 3]`
   :external:py:class:`numpy.ndarray`.
