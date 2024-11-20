fc := gfortran
fc_flags := -Ofast
mods := numbertypes.o customla_mod.o asvqrd_mod.o

all: test.out

test.out: test.f90 numbertypes.o customla_mod.o asvqrd_mod.o
	$(fc) $(fc_flags) test.f90 $(mods) -llapack -lblas -o test.out


numbertypes.o: numbertypes.f90
	$(fc) $(fc_flags) -c numbertypes.f90

customla_mod.o: customla_mod.f90 numbertypes.o
	$(fc) $(fc_flags) -c customla_mod.f90

asvqrd_mod.o: asvqrd_mod.f90 customla_mod.o numbertypes.o
	$(fc) $(fc_flags) -c asvqrd_mod.f90

clean:
	rm -rf *.o *.mod *.out