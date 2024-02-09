G4GOUPIl_DIR=$(shell python3 -m goupil --prefix)/interfaces/geant4

CFLAGS= -O2 \
        -I$(G4GOUPIl_DIR) \
        $(shell geant4-config --cflags)

LIBS= $(shell geant4-config --libs)

libgeometry.so: geometry.cpp
	$(CXX) $(CFLAGS) -shared -fPIC -o $@ $< $(G4GOUPIl_DIR)/G4Goupil.cc $(LIBS)

clean:
	rm -f libgeometry.so
