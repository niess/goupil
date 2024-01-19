.. _DistributionFunction:

`DistributionFunction`_
=======================

This class provides read-only views of tables that relate to the Cumulative
Distribution Function (CDF) of a collision process. The CDF values are tabulated
over an :math:`m \times n` grid of values, which depend on the energies of the
incoming and outgoing photons, :math:`\nu_i` and :math:`\nu_f`, respectively. To
map the support of :math:`\nu_f` onto :math:`(0,1)`, a logarithmic mapping is
used, as the support depends on the incoming energy :math:`\nu_i` (see the
:py:attr:`x <DistributionFunction.x>` attribute below).


Constructor
-----------

.. py:class:: DistributionFunction(*args, **kwargs)

.. warning::

   Distribution function views should not be instantiated directly, but only
   from a :doc:`material_record`. Direct instantiation will result in a
   :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: DistributionFunction.energies_in
   :type: numpy.ndarray

   The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for which
   the CDF was pre-computed.

.. py:attribute:: DistributionFunction.material
   :type: MaterialRecord

   The collision's target material.

.. py:attribute:: DistributionFunction.process
   :type: str

   A description of the interaction process.

.. py:attribute:: DistributionFunction.values
   :type: numpy.ndarray

   The :math:`n \times m` table of pre-computed CDF values.

.. py:attribute:: DistributionFunction.x
   :type: numpy.ndarray

   The :math:`n` mapped values for the energy :math:`\nu_f` of the outgoing
   photon, given as

   .. math::

      x = \frac{\ln(\nu_f/\nu_\text{min})}{\ln(\nu_\text{max}/\nu_\text{min})},

   where :math:`(\nu_\text{min}, \nu_\text{max})` is the DCS support, depending
   on :math:`\nu_i`. See the :py:meth:`energies_out` method below for the
   converse mapping.


Methods
-------

.. py:method:: DistributionFunction.__call__(energy_in, energy_out)

   Returns interpolated CDF value(s) for *energy_in* (:math:`\nu_i`) and
   *energy_out* (:math:`\nu_f`). The latter can be specified as a
   :external:py:class:`numpy.ndarray`.

.. py:method:: DistributionFunction.energies_out(index)

   Returns the :math:`m` values of the outgoing photon energy (:math:`\nu_f`)
   corresponding to the specified *index* for the incoming photon energy
   (:math:`\nu_i`).
