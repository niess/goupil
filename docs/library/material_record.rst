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

.. py:method:: MaterialRecord.compton_cdf(model=None, mode=None) -> CDF

   Returns a read-only view of the cumulative distribution function (CDF)
   describing a Compton collision.

.. py:method:: MaterialRecord.compton_cross_section(model=None, mode=None) -> CrossSection

   Returns a read-only view of the cross-section table for a given Compton
   process.

.. py:method:: MaterialRecord.compton_inverse_cdf(model=None, mode=None) -> InverseCDF

   Returns a read-only view of the inverse transform of the cumulative
   distribution function (CDF) describing a Compton collision.

.. py:method:: MaterialRecord.compton_weight(energy_in, energy_out, model=None, mode=None)

   This is a convenience function. It returns the Monte Carlo weight for the
   sampling of Compton collisions, for a given `model` and simulation `mode`.

.. py:method:: MaterialRecord.rayleigh_cross_section() -> CrossSection

   Returns a read-only view of the cross-section table for Rayleigh collisions.

.. py:method:: MaterialRecord.rayleigh_form_factor() -> FormFactor

   Returns a read-only view of the form factor describing Rayleigh collisions.


Table views
-----------

.. py:class:: CDF

   This class provides read-only views of tables that relate to the Cumulative
   Density Function (CDF) of a collision process. The CDF values are tabulated
   over an :math:`m \times n` grid of values, which depend on the energies of
   the incoming and outgoing photons, :math:`\nu_i` and :math:`\nu_f`,
   respectively. To map the support of :math:`\nu_f` onto :math:`(0,1)`, a
   logarithmic mapping is used, as it depends on the incoming energy
   :math:`\nu_i` (see the :py:attr:`x` attribute below).

   .. py:attribute:: energies_in
      :type: numpy.ndarray

      The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for
      which the CDF was pre-computed.

   .. py:attribute:: material
      :type: MaterialRecord

      The collision's target material.

   .. py:attribute:: process
      :type: str

      A description of the interaction process.

   .. py:attribute:: values
      :type: numpy.ndarray

      The :math:`n \times m` table of pre-computed CDF values.

   .. py:attribute:: x
      :type: numpy.ndarray

      The :math:`n` mapped values for the energy :math:`\nu_f` of the outgoing
      photon, given as

      .. math::

         x = \frac{\ln(\nu_f/\nu_\text{min})}{\ln(\nu_\text{max}/\nu_\text{min})},

      where :math:`(\nu_\text{min}, \nu_\text{max})` is the DCS support,
      depending on :math:`\nu_i`. See the :py:meth:`energies_out` method below
      for the converse mapping.

   .. py:method:: __call__(energy_in: float, energy_out: float | numpy.ndarray)

      Returns interpolated CDF value(s) for *energy_in* (:math:`\nu_i`) and
      *energy_out* (:math:`\nu_f`). The latter can be specified as a
      :external:py:class:`numpy.ndarray`.

   .. py:method:: energies_out(index: int)

      Returns the :math:`m` values of the outgoing photon energy (:math:`\nu_f`)
      corresponding to the specified *index* for the incoming photon energy
      (:math:`\nu_i`).


.. py:class:: CrossSection

   This class provides read-only views of tables that relate to the
   cross-section of a collision process. The cross-section values are tabulated
   over a grid of :math:`m` nodes, which depend on the energy :math:`\nu_i`
   of the incoming photon.

   .. py:attribute:: energies
      :type: numpy.ndarray

      The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for
      which the cross-section was pre-computed.

   .. py:attribute:: material
      :type: MaterialRecord

      The collision's target material.

   .. py:attribute:: process
      :type: str

      A description of the interaction process.

   .. py:attribute:: values
      :type: numpy.ndarray

      The :math:`m` pre-computed cross-section values.

   .. py:method:: __call__(energy: float | numpy.ndarray)

      Returns interpolated cross-section value(s) for the specified *energy*
      (:math:`\nu_i`) which can be either a float or a
      :external:py:class:`numpy.ndarray`.

.. py:class:: FormFactor

   This class provides read-only views of tables that relate to the form-factor
   of a collision process. The form-factor values are tabulated over a grid of
   :math:`m` nodes, which depend on the transfered momentum :math:`q`.

   .. py:attribute:: material
      :type: MaterialRecord

      The collision's target material.

   .. py:attribute:: momenta
      :type: numpy.ndarray

      The :math:`m` momenta values (:math:`q`), for which the form-factor was
      computed.

   .. py:attribute:: process
      :type: str

      A description of the interaction process.

   .. py:attribute:: values
      :type: numpy.ndarray

      The :math:`m` computed form-factor values.

   .. py:method:: __call__(momentum: float | numpy.ndarray)

      Returns interpolated form-factor value(s) for the specified *momentum*
      (:math:`q`) which can be either a float or a
      :external:py:class:`numpy.ndarray`.

.. py:class:: InverseCDF

   This class provides read-only views of tables that relate to the inverse of
   the Cumulative Density Function (CDF) of a collision process, i.e. energy
   values of the outgoing photon (:math:`\nu_f`). The inverse CDF values are
   tabulated over an :math:`m \times n` grid of values, which depend on the
   energy of the incoming photons (:math:`\nu_i`) and the CDF values (in
   :math:`(0,1)`).

   .. py:attribute:: cdf
      :type: numpy.ndarray

      The :math:`n` CDF values of the collision process, for
      which the outgoing photon energy (:math:`\nu_f`) was pre-computed.

   .. py:attribute:: energies
      :type: numpy.ndarray

      The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for
      which the inverse CDF was pre-computed.

   .. py:attribute:: material
      :type: MaterialRecord

      The collision's target material.

   .. py:attribute:: process
      :type: str

      A description of the interaction process.

   .. py:attribute:: values
      :type: numpy.ndarray

      The :math:`n \times m` table of pre-computed inverse CDF values, i.e.
      outgoing photon energies (:math:`\nu_f`).

   .. py:attribute:: weights
      :type: numpy.ndarray

      The :math:`n \times m` table of pre-computed sampling weights, in case of
      an inverse collision process, or :python:`None` otherwise.

   .. py:method:: __call__(energy: float, cdf: float | numpy.ndarray)

      Returns interpolated inverse CDF value(s) for the specified incoming
      *energy* (:math:`\nu_i`) and *cdf* value(s). The latter can be a
      :external:py:class:`numpy.ndarray`.
