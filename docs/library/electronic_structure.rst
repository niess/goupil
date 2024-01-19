.. _ElectronicStructure:

`ElectronicStructure`_
======================

This class represents the electronic structure of an atomic element, compound,
or mixture in terms of shells. It is a wrapper over a structured
:external:py:class:`numpy.ndarray` that contains shell data (see the
:py:attr:`shells <ElectronicStructure.shells>` method below).


Constructor
-----------

.. py:class:: ElectronicStructure(*args, **kwargs)

.. warning::

   Electronic structures cannot be instantiated directly, but only from an
   :doc:`atomic_element`, :doc:`material_definition` or :doc:`material_record`.
   Direct instantiation will result in a :external:py:class:`TypeError`.


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

   A structured :external:py:class:`numpy.ndarray` storing the properties of
   electronic shells, i.e. the shell binding *energy*, average *momentum* and
   *occupancy*. For example

   >>> H.electrons().shells["occupancy"]
   array([1.]
