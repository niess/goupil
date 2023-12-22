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

   .. note::

      If the requested parameters combination is not available, a
      :python:`NotImplementedError` is thrown.


Attributes
----------

.. py:attribute:: ComptonProcess.method
   :type: str

   This attribute determines the sampling method for Compton collisions.
   Available options are

   .. list-table::

      * - :python:`"Inverse Transform"`
        - Inverse transform sampling.

      * - :python:`"Rejection Sampling"`
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

   This attribute specifies the physical model describing Compton scattering.
   The possible values are

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

      When using the Klein-Nishina model, since the total cross-section is
      solved analytically, this parameter has no effect.


Methods
-------

.. py:method:: ComptonProcess.cross_section(energy, material, energy_min=None, energy_max=None)

   Computes the total cross-section for Compton scattering of a photon with a
   specified initial `energy` (in MeV) on an atom of a given `material`. The
   energy can be a :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray` of floats. The material must be an
   instance of :doc:`material_definition` or :doc:`material_record`. Optional
   bounds can be set on the energy of the outgoing photon using the `energy_min`
   and `energy_max` arguments (in MeV).

   Examples
   ^^^^^^^^

   >>> goupil.ComptonProcess().cross_section(1.0, H2O) # doctest: +SKIP

.. py:method:: ComptonProcess.dcs(energy_in, energy_out, material)

   Computes the differential cross-section (DCS) for Compton scattering of a
   photon on an atom of a given `material`. The input parameters are the ingoing
   energy (in MeV) and the outgoing energy (in MeV), which can be a single
   :external:py:class:`float` or a :external:py:class:`numpy.ndarray` of floats.
   The material must be provided as either a :doc:`material_definition` or
   :doc:`material_record` object.

   Examples
   ^^^^^^^^

   >>> goupil.ComptonProcess().dcs(1.0, 0.8, H2O) # doctest: +SKIP

.. py:method:: ComptonProcess.dcs_support(energy)

   Returns the support of the differential cross section (DCS) for a given
   photon `energy` in MeV. The energy can be specified as a
   :external:py:class:`float` or as a :external:py:class:`numpy.ndarray`. The
   output is a tuple containing the minimum and maximum bounds of the support.

   Examples
   ^^^^^^^^

   >>> goupil.ComptonProcess().dcs_support(1.0) # doctest: +SKIP

.. py:method:: ComptonProcess.sample(energy, material, rng=None)

   This function generates random Compton collisions. The input photon `energy`,
   in MeV, can be a :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray`. The target `material` must be an instance
   of :doc:`material_record`. The output is a tuple that contains the outgoing
   energy in MeV, the cosine of the scattering angle and the generation weight.
   It is also possible to provide a specific :doc:`random_stream` (`rng`) as an
   option.

   Examples
   ^^^^^^^^

   >>> goupil.ComptonProcess().sample(1.0, H20) # doctest: +SKIP
