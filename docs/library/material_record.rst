:class:`MaterialRecord`
=======================

.. _MaterialRecord:

----

This class stores tables of pre-calculated physical properties of a material
relevant to Monte Carlo transport, such as cross-sections. These data are
accessible as read-only :external:py:class:`numpy.ndarray` through table views
described below.


Constructor
-----------

.. py:class:: MaterialRecord(*args, **kwargs)

.. warning::

   Material records should not be instantiated directly, but only from a
   :doc:`material_registry`. Direct instantiation will result in a
   :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: Materialrecord.definition
   :type: MaterialDefinition

   The material definition from which the record was generated.

.. py:attribute:: Materialrecord.electrons
   :type: ElectronicStructure

   The electronic structure associated to this material.


Methods
-------

.. py:method:: MaterialRecord.absorption_cross_section() -> CrossSection

   Returns a read-only view of the cross-section table for absorption processes.

.. note::

   See :doc:`compton_process` for the meaning of `mode` and `model` parameters
   in Compton processes related methods below.

.. py:method:: MaterialRecord.compton_cdf(model=None, mode=None) -> DistributionFunction

   Returns a read-only view of the cumulative distribution function (CDF)
   describing a Compton collision.

.. py:method:: MaterialRecord.compton_cross_section(model=None, mode=None) -> CrossSection

   Returns a read-only view of the cross-section table for a given Compton
   process.

.. py:method:: MaterialRecord.compton_inverse_cdf(model=None, mode=None) -> InverseDistribution

   Returns a read-only view of the inverse transform of the cumulative
   distribution function (CDF) describing a Compton collision.

.. py:method:: MaterialRecord.compton_weight(energy_in, energy_out, model=None, mode=None)

   This is a convenience function. It returns the Monte Carlo weight for the
   sampling of Compton collisions, for a given `model` and simulation `mode`.

.. py:method:: MaterialRecord.rayleigh_cross_section() -> CrossSection

   Returns a read-only view of the cross-section table for Rayleigh collisions.

.. py:method:: MaterialRecord.rayleigh_form_factor() -> FormFactor

   Returns a read-only view of the form factor describing Rayleigh collisions.
