FC =mpif90
FFLAGS = -O3
#FFLAGS = -g -O0 
OBJECTS = ljnpt.o mcljnpt.o modljnpt.o ptljnpt.o obsvljnpt.o dsyev.o quenchf.o
LIBS = -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LIBS2 = -L/usr/local/lib -lstdc++
all: ljnpt.exe

ljnpt.exe:  $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(LIBS2)

ljnpt.o: ljnpt.f90
	$(FC) -c -o $@ $(FFLAGS) $^

mcljnpt.o: mcljnpt.f90
	$(FC) -c -o $@ $(FFLAGS) $^

modljnpt.o: modljnpt.f90
	$(FC) -c -o $@ $(FFLAGS) $^

ptljnpt.o : ptljnpt.f90
	$(FC) -c -o $@ $(FFLAGS) $^

obsvljnpt.o: obsvljnpt.f90
	$(FC) -c -o $@ $(FFLAGS) $^

dsyev.o: dsyev.f
	$(FC) -c -o $@ $(FFLAGS) $^

quenchf.o: quenchf.cpp
	g++ -c -o $@ $(FFLAGS) $^

mcljnpt.f90 : modljnpt.o
modljnpt.f90 : quenchmod.o

clean:
	$(RM) $(OBJECTS) ljnpt.exe

.PHONY: clean
