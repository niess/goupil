:class:`ElectronicStructure`
============================

.. _ElectronicStructure:

----

This class represents an electronic structure (XXX ellaborate)


Constructor
-----------

.. py:class:: ElectronicStructure(seed=None)


Attributes
----------

.. py:attribute:: ElectronicStructure.charge
   :type: float

   The effective charge number of the structure. For compounds, it is the sum of
   the charge numbers of its constituent parts, and for mixtures, it is the
   average charge number of its constituent parts. For instance

   >>> H2O.electrons().charge
   10.0

.. py:attribute:: ElectronicStructure.shells
   :type: numpy.ndarray

   A structured :external:py:class:`numpy.ndarray` providing the electronic
   shells binding energies, average momenta and occupancy numbers.

   >>> H.electrons().shells
   10.0
