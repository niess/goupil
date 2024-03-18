// Geant4 interface.
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"
#include "Randomize.hh"
// Goupil interface.
#include "G4Goupil.hh"


// ============================================================================
//
// Geant4 geometry implementation.
//
// ============================================================================

// Hard coded size parameters.
static const G4double WORLD_SIZE = 2.0 * CLHEP::km;
static const G4double DETECTOR_WIDTH = 20.0 * CLHEP::m;
static const G4double DETECTOR_HEIGHT = 10.0 * CLHEP::m;
static const G4double DETECTOR_OFFSET = 5.0 * CLHEP::cm;


struct DetectorConstruction: public G4VUserDetectorConstruction {
    G4VPhysicalVolume * Construct() {
        auto manager = G4NistManager::Instance();

        // World wolume, containing the atmosphere layer.
        G4LogicalVolume * world;
        {
            std::string name = "Atmosphere";
            auto solid = new G4Box(
                name,
                0.5 * WORLD_SIZE,
                0.5 * WORLD_SIZE,
                0.5 * WORLD_SIZE
            );
            auto material = manager->FindOrBuildMaterial("G4_AIR");
            world = new G4LogicalVolume(solid, material, name);
        }

        // Ground volume.
        {
            std::string name = "Soil";
            auto solid = new G4Box(
                name,
                0.5 * WORLD_SIZE,
                0.5 * WORLD_SIZE,
                0.25 * WORLD_SIZE
            );
            auto material = manager->FindOrBuildMaterial(
                "G4_CALCIUM_CARBONATE");
            auto volume = new G4LogicalVolume(solid, material, name);
            new G4PVPlacement(
                nullptr,
                G4ThreeVector(0.0, 0.0, -0.25 * WORLD_SIZE),
                volume,
                name,
                world,
                false,
                0
            );
        }

        // Collection volume.
        {
            std::string name = "Detector";
            auto solid = new G4Box(
                name,
                0.5 * DETECTOR_WIDTH,
                0.5 * DETECTOR_WIDTH,
                0.5 * DETECTOR_HEIGHT
            );
            auto material = manager->FindOrBuildMaterial("G4_AIR");
            auto volume = new G4LogicalVolume(solid, material, name);
            new G4PVPlacement(
                nullptr,
                G4ThreeVector(
                    0.0, 0.0, 0.5 * DETECTOR_HEIGHT + DETECTOR_OFFSET),
                volume,
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
            world->GetName(),
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


// ============================================================================
//
// Initialisation routines for Monte Carlo states.
//
// ============================================================================

static void InitialisePrng() {
    // Get a seed from /dev/urandom (assuming a Unix system).
    unsigned long seed = 0;
    FILE * fid = std::fopen("/dev/urandom", "rb");
    fread(&seed, sizeof(long), 1, fid);
    fclose(fid);

    // Initialize the Geant4 pseudo-random engine.
    G4Random::setTheEngine(new CLHEP::MTwistEngine);
    G4Random::setTheSeed(seed);
}


// LIbrary initialisation (using a constructor attribute).
__attribute__((constructor)) static void initialize() {
    InitialisePrng();
}


static double randomiseEnergy() {
    // Gamma emission lines (Pb-214 and Bi-214).
    // Index 0 is energy (in MeV) and index 1 is intensity (in percent).
    const int energy = 0, intensity = 1;
    const int n = 11;
    double source_spectrum[n][2] = {
        {0.242,  7.3},
        {0.295, 18.4},
        {0.352, 35.6},
        {0.609, 45.5},
        {0.768,  4.9},
        {0.934,  3.1},
        {1.120, 14.9},
        {1.238,  5.8},
        {1.378,  4.0},
        {1.764, 15.3},
        {2.204,  4.9},
    };

    // Compute the cumulative distribution.
    for (int i = 1; i < n; i++) {
        source_spectrum[i][intensity] += source_spectrum[i - 1][intensity];
    }

    // Randomise the emission line.
    const double r = G4UniformRand() * source_spectrum[n - 1][intensity];
    int i;
    for (i = 0; i < n - 1; i++) {
        if (r < source_spectrum[i][intensity]) break;
    }

    // Return the corresponding energy.
    return source_spectrum[i][energy];
}


static void randomiseForward(struct goupil_state * state) {
    // Randomise the emission line.
    state->energy = randomiseEnergy();

    // Randomise the emission direction, uniformly over the full solid angle.
    const double cos_theta = 2.0 * G4UniformRand() - 1.0;
    const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    const double phi = 2.0 * M_PI * G4UniformRand();
    state->direction.x = sin_theta * cos(phi);
    state->direction.y = sin_theta * sin(phi);
    state->direction.z = cos_theta;

    // Randomise the source position over the atmosphere, excluding the
    // detector volume.
    for (;;) {
        const G4double x = WORLD_SIZE * (G4UniformRand() - 0.5);
        const G4double y = WORLD_SIZE * (G4UniformRand() - 0.5);
        const G4double z = 0.5 * WORLD_SIZE * G4UniformRand();

        const G4double z0 = 0.5 * DETECTOR_HEIGHT + DETECTOR_OFFSET;
        if ((fabs(x) <= 0.5 * DETECTOR_WIDTH) &&
            (fabs(y) <= 0.5 * DETECTOR_WIDTH) &&
            (fabs(z - z0) <= 0.5 * DETECTOR_HEIGHT)) {
            // The tentative point lies inside the detector volume. Let us
            // generate another one.
            continue;
        }

        state->position.x = x / CLHEP::cm;
        state->position.y = y / CLHEP::cm;
        state->position.z = z / CLHEP::cm;
        break;
    }

    // Set the Monte Carlo weight.
    state->weight = 1.0;
}


extern "C" void initialise_states(
    size_t n, struct goupil_state * states, bool forward) {

    // Loop over Monte Carlo states assuming a *contiguous* memory layout.
    for (size_t i = 0; i < n; i++) {
        if (forward) {
            randomiseForward(states + i);
        }
    }
}
