:class:`AtomicElement`
======================

.. _AtomicElement:

----

This class represents an atomic element and provides access to its properties as
attributes. Atomic properties are immutable and predefined during the
initialization of :mod:`goupil`.


Constructor
-----------

.. py:class:: AtomicElement(arg: int | str)

   To create an instance of the desired element, use either its atomic number
   (as an :external:py:class:`int`) or its symbol (as a
   :external:py:class:`str`). For instance,

   >>> goupil.AtomicElement("H")
   H

   .. seealso::
      The :doc:`elements` function offers an alternative method for
      instantiating multiple atomic elements at once.


Attributes
----------

.. py:attribute:: AtomicElement.A
   :type: float

   The atomic element's mass number. For instance, hydrogen has a mass number of
   :python:`1.0087`.

.. py:attribute:: AtomicElement.name
   :type: str

   The complete name of the atomic element. For instance :python:`"hydrogen"`.

.. py:attribute:: AtomicElement.symbol
   :type: str

   The canonical symbol for the atomic element, i.e. a one or two-letter
   abbreviation. For instance, hydrogen is represented by the symbol
   :python:`"H"`.

.. py:attribute:: AtomicElement.Z
   :type: int

   The element's atomic number, represented by an :external:py:class:`int`. For
   instance, hydrogen has an atomic number of :python:`1`.
