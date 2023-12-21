:class:`ComptonProcess`
=======================

.. _compton_process:

----

This class provides an interface to a variety of implementations of the Compton
scattering process.

Constructor
-----------

.. py:class:: ComptonProcess(model: str=None, mode: str=None, method: str=None)

   The constructor arguments should match one of the attributes described below.


Attributes
----------

.. py:attribute:: ComptonProcess.method
   :type: str

   This attribute determines the sampling method for Compton collisions.
   Available options are

   .. list-table::

      * - :python:`"Inverse"`
        - Inverse transform sampling.

      * - :python:`"Rejection"`
        - Sampling by means of a rejection procedure.

   The inverse transform technique, the fastest, involves pre-computation and
   memory storage of differential cross sections that may become quite
   voluminous. Thus, by default, a rejection method is used which is
   approximately twice as slow as the inverse transform method.

.. py:attribute:: ComptonProcess.mode
   :type: str

   Compton collision simulation mode. The available options are as follows

   .. list-table::

      * - ``"Adjoint"``
        - Adjoint collision process.

      * - ``"Directe"``
        - Standard collision process.

      * - ``"Inverse"``
        - Reverse of the standard process.

      * - ``None``
        - Compton process is disabled.

   The default setting for reverse simulation is Adjoint, as it allows the use
   of a rejection method, unlike Inverse mode.


.. py:attribute:: ComptonProcess.model
   :type: str

   Specifies the physical model describing Compton scattering. The possible
   values are

   .. list-table::

      * - ``"ImpulseApproximation"``
        - The Impulse Approximation (IA) is used to account for the momentum
          distribution of the bound electrons, without any additional
          approximations.

      * - ``"KleinNishina"``
        - The electrons of the target atoms are assumed to be both free and at
          rest.

      * - ``"Penelope"``
        - The Penelope model [Baro95]_ is used.

      * - ``"ScatteringFunction"``
        - Effective model based on the Penelope scattering function [Baro95]_.

   The effective model, based on Penelope's Compton scattering function
   [Baro95]_, is used by default.
