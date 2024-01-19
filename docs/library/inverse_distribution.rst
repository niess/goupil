.. _InverseDistribution:

`InverseDistribution`_
======================

----

This class provides read-only views of tables that relate to the inverse of the
Cumulative Density Function (CDF) of a collision process, i.e. energy values of
the outgoing photon (:math:`\nu_f`). The inverse CDF values are tabulated over
an :math:`m \times n` grid of values, which depend on the energy of the incoming
photons (:math:`\nu_i`) and the CDF values (in :math:`(0,1)`).


Constructor
-----------

.. py:class:: InverseDistribution(*args, **kwargs)

.. warning::

   Inverse CDF views should not be instantiated directly, but only from a
   :doc:`material_record`. Direct instantiation will result in a
   :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: InverseDistribution.cdf
   :type: numpy.ndarray

   The :math:`n` CDF values of the collision process, for which the outgoing
   photon energy (:math:`\nu_f`) was pre-computed.

.. py:attribute:: InverseDistribution.energies
   :type: numpy.ndarray

   The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for which
   the inverse CDF was pre-computed.

.. py:attribute:: InverseDistribution.material
   :type: MaterialRecord

   The collision's target material.

.. py:attribute:: InverseDistribution.process
   :type: str

   A description of the interaction process.

.. py:attribute:: InverseDistribution.values
   :type: numpy.ndarray

   The :math:`n \times m` table of pre-computed inverse CDF values, i.e.
   outgoing photon energies (:math:`\nu_f`).

.. py:attribute:: InverseDistribution.weights
   :type: numpy.ndarray

   The :math:`n \times m` table of pre-computed sampling weights, in case of an
   inverse collision process, or :python:`None` otherwise.


Methods
-------

.. py:method:: InverseDistribution.__call__(energy, cdf)

   Returns interpolated inverse CDF value(s) for the specified incoming *energy*
   (:math:`\nu_i`) and *cdf* value(s). The latter can be a
   :external:py:class:`numpy.ndarray`.
