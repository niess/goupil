#ifndef goupil_h
#define goupil_h
#ifdef __cplusplus
extern "C" {
#endif

/* C standard library. */
#include <stddef.h>


/* ============================================================================
 *  External geometry interface.
 * ============================================================================
 */
struct goupil_interface {
    struct goupil_geometry_definition * (*new_geometry_definition)(void);

    struct goupil_geometry_tracer * (*new_geometry_tracer)(
        struct goupil_geometry_definition * definition
    );
};

struct goupil_interface goupil_initialise(void);


/* ============================================================================
 *  Float interface.
 * ============================================================================
 */
typedef double goupil_float_t;

struct goupil_float3 {
    goupil_float_t x;
    goupil_float_t y;
    goupil_float_t z;
};


/* ============================================================================
 *  Monte Carlo state interface.
 * ============================================================================
 */
struct goupil_state {
    goupil_float_t energy;
    struct goupil_float3 position;
    struct goupil_float3 direction;
    goupil_float_t length;
    goupil_float_t weight;
};


/* ============================================================================
 *  Geometry definition interface.
 * ============================================================================
 */
struct goupil_geometry_definition {
    /* Destroys the geometry. */
    void (*destroy)(struct goupil_geometry_definition * self);

    /* Returns the definition of an indexed material. */
    const struct goupil_material_definition * (*get_material)(
        const struct goupil_geometry_definition * self,
        size_t index
    );

    /* Returns data relative to a specific geometry sector. */
    const struct goupil_geometry_sector (*get_sector)(
        const struct goupil_geometry_definition * self,
        size_t index
    );

    /* Returns the total number of materials for this geometry. */
    size_t (*materials_len)(const struct goupil_geometry_definition * self);

    /* Returns the total number of sectors composing this geometry. */
    size_t (*sectors_len)(const struct goupil_geometry_definition * self);
};


/* ============================================================================
 *  Geometry tracer interface.
 * ============================================================================
 */
struct goupil_geometry_tracer {
    const struct goupil_geometry_definition * definition;

    void (*destroy)(struct goupil_geometry_tracer * self);

    struct goupil_float3 (*position)(
        const struct goupil_geometry_tracer * self
    );

    void (*reset)(
        struct goupil_geometry_tracer * self,
        struct goupil_float3 position,
        struct goupil_float3 direction
    );

    size_t (*sector)(const struct goupil_geometry_tracer * self);

    goupil_float_t (*trace)(
        struct goupil_geometry_tracer * self,
        goupil_float_t physical_length
    );

    void (*update)(
        struct goupil_geometry_tracer * self,
        goupil_float_t length,
        struct goupil_float3 direction
    );
};


/* ============================================================================
 *  Geometry sector interface.
 * ============================================================================
 */
struct goupil_geometry_sector {
    /* Index of constitutive material. */
    size_t material;

    /* Bulk density of this geometry sector. */
    goupil_float_t density;

    /* Brief definition of this geometry sector. */
    const char * definition;
};


/* ============================================================================
 *  Material definition interface.
 * ============================================================================
 */
struct goupil_material_definition {
    /* Returns the number of atomic elements composing this material. */
    size_t (*composition_len)(const struct goupil_material_definition * self);

    /* Returns data relative to a specific atomic element. */
    const struct goupil_weighted_element (*get_composition)(
        const struct goupil_material_definition * self,
        size_t index
    );

    /* Returns the material name. */
    const char * (*name)(const struct goupil_material_definition * self);
};

struct goupil_weighted_element {
    /* Atomic number of this element. */
    int Z;

    /* Molar weight of this element. */
    goupil_float_t weight;
};


#ifdef __cplusplus
}
#endif
#endif
