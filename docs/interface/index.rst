C interface
===========

.. _c_interface:

This section provides information on how to interface an existing geometry
engine with Goupil, at the C level. To adapt a specific geometry engine to
Goupil, one needs to create a concrete implementation of Goupil's C interface,
which is thus called an *adapter* herein. The `description`_ section below
offers insight into how to implement an adapter, while the `functions`_ and
`types`_ sections provide a technical overview.

.. note::

   Goupil comes with a source distribution of a `Geant4
   <https://geant4.web.cern.ch/>`_ adapter, which also serves as a detailed
   example of an adapter implementation. For instructions on how to use this
   adapter, see the :doc:`../description` section.

.. warning::

   Prior knowledge of Goupil's geometry model is assumed in this section. To aid
   comprehension, it is advisable to first become acquainted with
   :doc:`../library/external_geometry` objects through the Python interface.


Description
-----------

The adapter's main purpose is to provide Goupil with a Monte Carlo geometry that
has ray tracing capabilities. This geometry should be packaged in a shared
library, requiring an entry point at runtime. The entry point is defined as the
:c:func:`goupil_initialise <goupil_initialise>` function.

.. note::

   The :c:func:`goupil_initialise <goupil_initialise>` function must be
   implemented exactly as written to ensure proper location by Goupil when
   loading the shared library.

The adapter must implement two crucial objects to interact with the Monte Carlo
geometry. The first object, :doc:`geometry_definition`, is used by Goupil to
build a replica of the Monte Carlo geometry using its internal representation.
The adapter must ensure that the associated Monte Carlo geometry remains
unchanged and active throughout the lifetime of the definition object.
Conversely, when the definition object is destroyed, the Monte Carlo geometry is
no longer used.

The second object, :doc:`geometry_tracer`, is used by Goupil to navigate through
the geometry during Monte Carlo steps. The adapter must conform to a
step-by-step algorithm, starting with a :c:func:`reset
<goupil_geometry_tracer.reset>` call, and then chaining sequences of
:c:func:`trace <goupil_geometry_tracer.trace>` and :c:func:`update
<goupil_geometry_tracer.update>` calls.

.. note::

   Unlike the :doc:`geometry_definition` object, tracer objects are volatile.
   This means that multiple tracer objects may be created and destroyed during
   the lifetime of a definition object. However, the adapter's implementation
   must ensure that the associated geometry object remains unchanged and alive
   throughout the lifetime of the tracers.

Goupil's C interface defines additional types that act as containers or
interfaces to containers. Please refer to the `types`_ section below for a
complete list.


Functions
---------

.. toctree::
   :maxdepth: 1

   initialise


Types
-----

.. toctree::
   :maxdepth: 1

   float3
   float_t
   geometry_definition
   geometry_sector
   geometry_tracer
   interface
   material_definition
   weighted_element
