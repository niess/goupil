:class:`ComptonProcess`
=======================

.. _compton_process:

----

This class provides an interface to a variety of implementations of the Compton
scattering process.

Constructor
-----------

.. py:class:: ComptonProcess(**kwargs)

   The constructor argument(s) should match one of the attributes described
   below. For instance

   >>> goupil.ComptonProcess(model="Penelope", precision=10.0)
   ComptonProcess(model="Penelope", precision=10)


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

      * - :python:`"Adjoint"`
        - Adjoint collision process.

      * - :python:`"Directe"`
        - Standard collision process.

      * - :python:`"Inverse"`
        - Reverse of the standard process.

      * - :python:`None`
        - Compton process is disabled.

   The default setting for reverse simulation is Adjoint, as it allows the use
   of a rejection method, unlike Inverse mode.


.. py:attribute:: ComptonProcess.model
   :type: str

   Specifies the physical model describing Compton scattering. The possible
   values are

   .. list-table::
      :widths: 26 59

      * - :python:`"Klein-Nishina"`
        - The electrons of the target atoms are assumed to be both free and at
          rest.

      * - :python:`"Penelope"`
        - The Penelope model [Baro95]_ is used, where the Impulse Approximation
          (IA) accounts for the momentum distribution of the bound electrons.

      * - :python:`"Scattering Function"`
        - Effective model based on the Penelope scattering function [Baro95]_.

   The effective model, based on Penelope's Compton scattering function
   [Baro95]_, is used by default.


.. py:attribute:: ComptonProcess.precision
   :type: float

   The numeric precision for cross-section computations relative to 1, which is
   the default value.

   .. note::

      When using the Klein-Nishina model, since the total cross-section can be
      solved analytically, this parameter has no effect.
