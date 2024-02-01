.. _RayleighProcess:

`RayleighProcess`_
==================

This class provides access to the implementation of the Rayleigh scattering
process, as detailed in [Baro95]_.

.. warning::

   All methods described below are defined at the class level.


Constructor
-----------

.. py:class:: RayleighProcess

   .. note::

      This class serves as a namespace and is not intended for instantiation.
      Any attempt to do so will result in a :external:py:class:`TypeError`.


Methods
-------

.. py:method:: RayleighProcess.cross_section(energy, material)

   Computes the total cross-section for Rayleigh scattering of a photon with a
   specified *energy* (in MeV) on an atom of a given *material*. The *energy*
   can be a :external:py:class:`float` or a :external:py:class:`numpy.ndarray`
   of floats. The *material* must be consistent with a
   :doc:`material_definition`. For instance,

   >>> goupil.RayleighProcess.cross_section(0.1, "H2O")
   1.600238e-25

.. py:method:: RayleighProcess.dcs(energy, cos_theta, material)

   Computes the differential cross-section (DCS) for Rayleigh scattering of a
   photon on an atom of a given *material*. The input parameters are the photon
   *energy* (in MeV) and the cosine of the scattering angle (*cos_theta*), which
   can be a single :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray` of floats. The *material* must be
   consistent with a :doc:`material_definition`.

.. py:method:: RayleighProcess.sample(energy, material, rng=None)

   This function generates random Rayleigh collisions. The photon *energy*, in
   MeV, can be a :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray`. The target *material* must be consistent
   with a :doc:`material_definition`. The output is the cosine of the scattering
   angle(s). It is also possible to provide a specific :doc:`random_stream`
   (`rng`) as an option.
