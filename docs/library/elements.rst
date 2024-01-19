.. _elements:

`elements`_
===========

This is a utility function that is used to instantiate multiple
:doc:`atomic_element` elements at the same time.

Syntax
------

.. py:function:: elements(*args: int | str)

   The input for this function consists of either the atomic numbers
   (:external:py:class:`int`) or symbols (:external:py:class:`str`) of the
   desired elements. It is also possible to specify a sequence of symbols
   directly as a single character string. Please refer to the `examples`_ below
   for further clarification.


Examples
--------

.. _examples:

Instantiation using atomic numbers.

>>> goupil.elements(1, 8)
(H, O)

Instanciation using atomic symbols.

>>> goupil.elements("H", "O")
(H, O)

Alternatively, elements may be instantiated from a literal sequence of atomic
symbols.

>>> goupil.elements("H, O")
(H, O)
