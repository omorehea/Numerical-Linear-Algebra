FC = gfortran

FFLAGS = -Wall -Wextra -Wimplicit-interface -Wno-surprising -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -fdefault-real-8 -fdefault-double-8

OBJECTS = utility.o LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

LinAl: $(OBJECTS)
	$(FC) -o $@ $^ $(FFLAGS)

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	$(FC) -c $< $(FFLAGS)

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<
.PHONY: clean

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl  LinAl.exe *.o *.mod *~
