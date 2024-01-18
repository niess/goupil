goupil_initialise_
==================

.. _goupil_initialise:

----

Initialises a C interface to Goupil.


Prototype
---------

.. c:function:: struct goupil_interface goupil_initialise(void)


Description
-----------

The main entry point in a shared library that exposes a Monte Carlo geometry to
Goupil is this function. Upon completion, it should return a ready-to-use
:c:struct:`goupil_interface` structure.
