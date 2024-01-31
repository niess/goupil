.. _goupil_state:

`goupil_state`_
===============

This structure contains a Monte Carlo state representing a transported photon.

.. c:struct:: goupil_state

   .. c:var:: goupil_float_t energy

      The photon energy, in MeV.

   .. c:var:: goupil_float3 position

      The photon position, in cm.

   .. c:var:: goupil_float3 direction

      The photon momentum's direction unit vector.

   .. c:var:: goupil_float_t length

      The photon total path length.

   .. c:var:: goupil_float_t weight

      The photon Monte Carlo weight.
