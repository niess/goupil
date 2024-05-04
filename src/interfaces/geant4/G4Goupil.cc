// Geant4 interface.
#include "G4Navigator.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
// Goupil C interface.
#include "goupil.h"
// Goupil Geant4 interface.
#include "G4Goupil.hh"
// C++ standard library.
#include <unordered_map>

// Entry point for Goupil.
#ifndef G4GOUPIL_INITIALISE
#define G4GOUPIL_INITIALISE goupil_initialise
#endif


// ============================================================================
//
// Local interface, bridging Geant4 and Goupil.
//
// ============================================================================

namespace G4Goupil {
    struct MaterialDefinition;

    class GeometryDefinition: public goupil_geometry_definition {
        public:
            GeometryDefinition(const G4VPhysicalVolume * world);
            ~GeometryDefinition();

            const struct goupil_material_definition *
                GetMaterial(size_t index) const;

            size_t GetMaterialIndex(const G4Material * material) const;

            const struct goupil_geometry_sector GetSector(size_t index) const;

            size_t GetSectorIndex(const G4VPhysicalVolume * volume) const;

            size_t MaterialsLen() const;

            size_t SectorsLen() const;

            const G4VPhysicalVolume * GetWorld() const;

        private:
            std::vector<const MaterialDefinition *> materials;
            std::vector<const G4VPhysicalVolume *> volumes;
            std::unordered_map<const G4Material *, size_t> materialsIndices;
            std::unordered_map<const G4VPhysicalVolume *, size_t>
                sectorsIndices;
    };

    struct GeometryTracer: public goupil_geometry_tracer {
        GeometryTracer(const GeometryDefinition * geometry);
        ~GeometryTracer();

        // Internal data.
        G4ThreeVector currentDirection;
        size_t currentIndex;
        G4ThreeVector currentPosition;
        goupil_float_t stepLength;
        goupil_float_t stepSafety;

        G4TouchableHistory * history;
        G4Navigator navigator;
    };

    struct MaterialDefinition: public goupil_material_definition {
        MaterialDefinition(const G4Material * material);
        ~MaterialDefinition();

        const G4Material * g4Material;
    };
}

// ============================================================================
//
// Implementation of Goupil C interface.
//
// ============================================================================

static struct goupil_geometry_definition * new_geometry_definition(void) {
    auto && topVolume = G4Goupil::NewGeometry();
    return new G4Goupil::GeometryDefinition(topVolume);
}


static struct goupil_geometry_tracer * new_geometry_tracer(
    struct goupil_geometry_definition * definition
){
    auto geometry = (G4Goupil::GeometryDefinition *)definition;
    return new G4Goupil::GeometryTracer(geometry);
}


extern "C" struct goupil_interface G4GOUPIL_INITIALISE (void) {
    struct goupil_interface interface;
    interface.new_geometry_definition = &new_geometry_definition;
    interface.new_geometry_tracer = &new_geometry_tracer;
    return interface;
}


// ============================================================================
//
// Implementation of geometry definition.
//
// ============================================================================

static void geometry_destroy(struct goupil_geometry_definition * self) {
    auto geometry = (G4Goupil::GeometryDefinition *)self;
    G4Goupil::DropGeometry(geometry->GetWorld());
    delete geometry;
}


static const struct goupil_material_definition * geometry_get_material(
    const struct goupil_geometry_definition * self,
    size_t index
){
    auto geometry = (G4Goupil::GeometryDefinition *)self;
    return geometry->GetMaterial(index);
}


static const struct goupil_geometry_sector geometry_get_sector(
    const struct goupil_geometry_definition * self,
    size_t index
){
    auto geometry = (G4Goupil::GeometryDefinition *)self;
    return geometry->GetSector(index);
}


static size_t geometry_materials_len(
    const struct goupil_geometry_definition * self
){
    auto geometry = (G4Goupil::GeometryDefinition *)self;
    return geometry->MaterialsLen();
}


static size_t geometry_sectors_len(
const struct goupil_geometry_definition * self)
{
    auto geometry = (G4Goupil::GeometryDefinition *)self;
    return geometry->SectorsLen();
}


static void append(
    std::vector<const G4Goupil::MaterialDefinition *> &materials,
    std::vector<const G4VPhysicalVolume *> &volumes,
    std::unordered_map<const G4Material *, size_t> &materialsIndices,
    std::unordered_map<const G4VPhysicalVolume *, size_t> &sectorsIndices,
    const G4VPhysicalVolume * current
){
    if (sectorsIndices.count(current) == 0) {
        size_t n = volumes.size();
        sectorsIndices.insert({current, n});
        volumes.push_back(current);

        auto && material = current->GetLogicalVolume()->GetMaterial();
        if (materialsIndices.count(material) == 0) {
            size_t m = materials.size();
            materialsIndices.insert({material, m});
            materials.push_back(new G4Goupil::MaterialDefinition(material));
        }
    }
    auto && logical = current->GetLogicalVolume();
    G4int n = logical->GetNoDaughters();
    for (G4int i = 0; i < n; i++) {
        auto && volume = logical->GetDaughter(i);
        append(materials, volumes, materialsIndices, sectorsIndices, volume);
    }
}


G4Goupil::GeometryDefinition::GeometryDefinition(
    const G4VPhysicalVolume * world)
{
    // Set interface.
    this->destroy = &geometry_destroy;
    this->get_material = &geometry_get_material;
    this->get_sector = &geometry_get_sector;
    this->materials_len = &geometry_materials_len;
    this->sectors_len = &geometry_sectors_len;

    // Scan volumes hierarchy.
    append(
        this->materials,
        this->volumes,
        this->materialsIndices,
        this->sectorsIndices,
        world
    );
}


G4Goupil::GeometryDefinition::~GeometryDefinition() {
    for (auto && material: this->materials) {
        delete material;
    }
}


const struct goupil_material_definition *
    G4Goupil::GeometryDefinition::GetMaterial(size_t index) const {
    if (index < this->materials.size()) {
        return this->materials[index];
    } else {
        return nullptr;
    }
}


size_t G4Goupil::GeometryDefinition::GetMaterialIndex(
    const G4Material * material) const
{
    try {
        return this->materialsIndices.at(material);
    } catch (...) {
        return this->materialsIndices.size();
    }
}


const struct goupil_geometry_sector
    G4Goupil::GeometryDefinition::GetSector(size_t index) const
{
    if (index < this->volumes.size()) {
        auto && volume = this->volumes[index];
        auto && material = volume->GetLogicalVolume()->GetMaterial();
        size_t material_index = this->GetMaterialIndex(material);
        return {
            material_index,
            material->GetDensity() / CLHEP::g * CLHEP::cm3,
            volume->GetName().c_str()
        };
    } else {
        return {
            this->materials.size(),
            0.0,
            nullptr
        };
    }
}


size_t G4Goupil::GeometryDefinition::GetSectorIndex(
    const G4VPhysicalVolume * volume) const
{
    try {
        return this->sectorsIndices.at(volume);
    } catch (...) {
        return this->sectorsIndices.size();
    }
}


size_t G4Goupil::GeometryDefinition::MaterialsLen() const {
    return this->materials.size();
}


size_t G4Goupil::GeometryDefinition::SectorsLen() const {
    return this->volumes.size();
}


const G4VPhysicalVolume * G4Goupil::GeometryDefinition::GetWorld() const {
    return (this->volumes.size() > 0) ?
        this->volumes[0] :
        nullptr;
}


// ============================================================================
//
// Implementation of geometry tracer.
//
// ============================================================================

static void tracer_destroy(struct goupil_geometry_tracer * self) {
    auto tracer = (G4Goupil::GeometryTracer *)self;
    delete tracer;
}


static struct goupil_float3 tracer_position(
    const struct goupil_geometry_tracer * self
){
    auto tracer = (G4Goupil::GeometryTracer *)self;
    auto && r = tracer->currentPosition;
    struct goupil_float3 position_ = {
        (goupil_float_t)(r[0] / CLHEP::cm),
        (goupil_float_t)(r[1] / CLHEP::cm),
        (goupil_float_t)(r[2] / CLHEP::cm)
    };
    return position_;
}


static void tracer_reset(
    struct goupil_geometry_tracer * self,
    struct goupil_float3 position,
    struct goupil_float3 direction
){
    auto tracer = (G4Goupil::GeometryTracer *)self;
    auto geometry = (const G4Goupil::GeometryDefinition *)self->definition;

    // Reset Geant4 navigation.
    tracer->currentPosition = G4ThreeVector(
        position.x * CLHEP::cm,
        position.y * CLHEP::cm,
        position.z * CLHEP::cm
    );

    tracer->currentDirection = G4ThreeVector(
        direction.x,
        direction.y,
        direction.z
    );

    tracer->navigator.ResetStackAndState();
    tracer->navigator.LocateGlobalPointAndUpdateTouchable(
        tracer->currentPosition,
        tracer->currentDirection,
        tracer->history,
        false // Do not use history.
    );

    // Reset internal state.
    tracer->currentIndex = geometry->GetSectorIndex(
        tracer->history->GetVolume()
    );
    tracer->stepLength = 0.0;
    tracer->stepSafety = 0.0;
}


static size_t tracer_sector(const struct goupil_geometry_tracer * self){
    auto tracer = (G4Goupil::GeometryTracer *)self;
    return tracer->currentIndex;
}


static goupil_float_t tracer_trace(
    struct goupil_geometry_tracer * self,
    goupil_float_t physical_length
){
    auto tracer = (G4Goupil::GeometryTracer *)self;

    G4double safety = 0.0;
    G4double s = tracer->navigator.ComputeStep(
        tracer->currentPosition,
        tracer->currentDirection,
        physical_length * CLHEP::cm,
        safety
    );
    goupil_float_t step = (goupil_float_t)(s / CLHEP::cm);
    tracer->stepLength = step;
    tracer->stepSafety = (goupil_float_t)(safety / CLHEP::cm);

    return (step < physical_length) ? step : physical_length;
}


static void tracer_update(
    struct goupil_geometry_tracer * self,
    goupil_float_t length,
    struct goupil_float3 direction
){
    auto tracer = (G4Goupil::GeometryTracer *)self;

    tracer->currentPosition += (length * CLHEP::cm) * tracer->currentDirection;

    if (length < tracer->stepSafety) {
        tracer->navigator.LocateGlobalPointWithinVolume(
            tracer->currentPosition
        );
    } else {
        if (length >= tracer->stepLength) {
            tracer->navigator.SetGeometricallyLimitedStep();
        }
        tracer->navigator.LocateGlobalPointAndUpdateTouchable(
            tracer->currentPosition,
            tracer->currentDirection,
            tracer->history
        );
        auto geometry = (const G4Goupil::GeometryDefinition *)self->definition;
        tracer->currentIndex = geometry->GetSectorIndex(
            tracer->history->GetVolume()
        );
    }

    tracer->currentDirection = G4ThreeVector(
        direction.x,
        direction.y,
        direction.z
    );
}


G4Goupil::GeometryTracer::GeometryTracer(
    const G4Goupil::GeometryDefinition * definition_
){
    // Initialise Geant4 navigator.
    this->navigator.SetWorldVolume(
        (G4VPhysicalVolume *) definition_->GetWorld());
    this->history = this->navigator.CreateTouchableHistory();

    // Initialise internal data.
    this->currentDirection = G4ThreeVector(0.0, 0.0, 1.0);
    this->currentIndex = 0;
    this->currentPosition = G4ThreeVector(0.0, 0.0, 0.0);
    this->stepLength = 0.0;
    this->stepSafety = 0.0;

    // Set C interface.
    this->definition = (const goupil_geometry_definition *)definition_;

    this->destroy = &tracer_destroy;
    this->position = &tracer_position;
    this->reset = &tracer_reset;
    this->sector = &tracer_sector;
    this->trace = &tracer_trace;
    this->update = &tracer_update;
}


G4Goupil::GeometryTracer::~GeometryTracer() {
    delete this->history;
}


// ============================================================================
//
// Implementation of material definition.
//
// ============================================================================

static size_t material_composition_len(
    const struct goupil_material_definition * self
){
    auto material = (G4Goupil::MaterialDefinition *)self;
    return material->g4Material->GetNumberOfElements();
}


static const struct goupil_weighted_element material_get_composition(
    const struct goupil_material_definition * self,
    size_t index
){
    auto material = ((G4Goupil::MaterialDefinition *)self)->g4Material;
    if (index < material->GetNumberOfElements()) {
        int Z = int(material->GetElement(index)->GetZ());
        goupil_float_t weight = goupil_float_t(
            material->GetVecNbOfAtomsPerVolume()[index] /
            material->GetTotNbOfAtomsPerVolume()
        );
        return { Z, weight };
    } else {
        return { 0, 0.0 };
    }
}


static const char * material_name(
    const struct goupil_material_definition * self
){
    auto material = (G4Goupil::MaterialDefinition *)self;
    return material->g4Material->GetName().c_str();
}


G4Goupil::MaterialDefinition::MaterialDefinition(
    const G4Material * material
):
    g4Material(material)
{
    // Set interface.
    this->composition_len = &material_composition_len;
    this->get_composition = &material_get_composition;
    this->name = &material_name;
}

G4Goupil::MaterialDefinition::~MaterialDefinition() {}
