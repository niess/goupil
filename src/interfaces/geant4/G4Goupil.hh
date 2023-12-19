#ifndef G4Goupil_hh
#define G4Goupil_hh

class G4VPhysicalVolume;

namespace G4Goupil {
    const G4VPhysicalVolume * NewGeometry();
    void DropGeometry(const G4VPhysicalVolume * volume);
}

#endif
