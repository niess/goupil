// Geant4 interface.
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"
// Goupil interface.
#include "G4Goupil.hh"


// ============================================================================
//
// Geant4 geometry implementation.
//
// ============================================================================

static const G4double SIZE = 2.0 * CLHEP::km;
static const G4double INNER_SIZE = 20.0 * CLHEP::m;

struct DetectorConstruction: public G4VUserDetectorConstruction {
    G4VPhysicalVolume * Construct() {
        auto manager = G4NistManager::Instance();

        // World wolume, containing the atmosphere layer.
        G4LogicalVolume * world;
        {
            std::string name = "Atmosphere";
            auto solid = new G4Box(name, 0.5 * SIZE, 0.5 * SIZE, 0.5 * SIZE);
            auto material = manager->FindOrBuildMaterial("G4_AIR");
            world = new G4LogicalVolume(solid, material, name);
        }

        // Ground volume.
        {
            std::string name = "Ground";
            auto solid = new G4Box(name, 0.5 * SIZE, 0.5 * SIZE, 0.25 * SIZE);
            auto material = manager->FindOrBuildMaterial(
                "G4_CARBONATE_CALCIUM");
            auto ground = new G4LogicalVolume(solid, material, name);
            G4PVPlacement(
                nullptr,
                G4ThreeVector(0.0, 0.0, -0.5 * SIZE),
                ground,
                name,
                world,
                false,
                0
            );
        }

        // Inner volume.
        {
            std::string name = "Ground";
            auto solid = new G4Box(name, 0.5 * INNER_SIZE, 0.5 * INNER_SIZE,
                0.25 * INNER_SIZE);
            auto material = manager->FindOrBuildMaterial("G4_AIR");
            auto inner = new G4LogicalVolume(solid, material, name);
            G4PVPlacement(
                nullptr,
                G4ThreeVector(0.0, 0.0, 0.25 * INNER_SIZE),
                inner,
                name,
                world,
                false,
                0
            );
        }

        return new G4PVPlacement(
            nullptr,
            G4ThreeVector(0.0, 0.0, 0.0),
            world,
            "World",
            nullptr,
            false,
            0
        );
    }
};


// ============================================================================
//
// Goupil hooks.
//
// ============================================================================

const G4VPhysicalVolume * G4Goupil::NewGeometry() {
    // Build the geometry and return the top "World" volume.
    return DetectorConstruction().Construct();
}

void G4Goupil::DropGeometry(const G4VPhysicalVolume * volume) {
    // Delete any sub-volume(s).
    auto && logical = volume->GetLogicalVolume();
    while (logical->GetNoDaughters()) {
        auto daughter = logical->GetDaughter(0);
        logical->RemoveDaughter(daughter);
        G4Goupil::DropGeometry(daughter);
    }
    // Delete this volume.
    delete logical->GetSolid();
    delete logical;
    delete volume;
}
