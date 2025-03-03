default: elpot

SHTOOLS = "*directory with SHTOOLS library*"

FC = ifx
FFLAGS =
FFLAGS2=-m64 -fpp -free -O3 -Tf 
LNK=-I$(SHTOOLS)/include -L$(SHTOOLS)/lib -lSHTOOLS -L/usr/local -lfftw3 -lm -qmkl

obj/Parameters.o: Parameters.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c Parameters.f90      -o obj/Parameters.o      -module obj   
obj/Potential.o: Potential.f90
	$(FC) $(FFLAGS) -c Potential.f90       -o obj/Potential.o       -module obj
obj/MathLibrary.o: MathLibrary.f90
	$(FC) $(FFLAGS) -c MathLibrary.f90     -o obj/MathLibrary.o     -module obj
obj/Miller.o: Miller.f90
	$(FC) $(FFLAGS) -c Miller.f90          -o obj/Miller.o          -module obj
obj/Qshep.o: Qshep.f90
	$(FC) $(FFLAGS) -c Qshep.f90           -o obj/Qshep.o           -module obj
obj/DMol_files.o: DMol_files.f90
	$(FC) $(FFLAGS) -c DMol_files.f90      -o obj/DMol_files.o      -module obj
obj/QE_files.o: QE_files.f90
	$(FC) $(FFLAGS) -c QE_files.f90        -o obj/QE_files.o        -module obj
obj/VASP_files.o: VASP_files.f90
	$(FC) $(FFLAGS) -c VASP_files.f90      -o obj/VASP_files.o      -module obj
obj/TecPlot_module.o: TecPlot_module.f90
	$(FC) $(FFLAGS) -c TecPlot_module.f90  -o obj/TecPlot_module.o  -module obj


OBJS =    obj/Parameters.o     \
          obj/Potential.o      \
          obj/Miller.o         \
          obj/MathLibrary.o    \
          obj/Qshep.o  \
          obj/DMol_files.o     \
          obj/VASP_files.o     \
          obj/QE_files.o       \
          obj/TecPlot_module.o 


elpot: $(OBJS)
	$(FC) $(FFLAGS2) SH_analysis.f90 main.f90 $(OBJS) $(LNK) -o ~/bin/elpot -module obj


clean:
	rm -rf ~/bin/elpot obj/*.o obj/*.mod


