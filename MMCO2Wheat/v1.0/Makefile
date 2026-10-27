# Program name (name of the executable file)
PROG = mmco2crop

# Object files (TYPES.o must be first as others depend on it)
OBJS = TYPES.o SHARED.o csv_file.o radiation_module.o soil_heat_module.o energy_balance_module.o CANOPY_MODEL.o growth_module.o mineral_module.o soiln_module.o phenology_module.o main_MMCO2Crop.o count_lines.o th.o meteo.o solar.o hydro.o assim.o

# Compiler and Flags
FC     = gfortran
FFLAGS = -I/usr/local/include 
FLINK  = -L/usr/local/lib
LINKER = $(FC) -o
switch = -O3 -Wall

model: $(PROG)

# Dependencies
TYPES.o: TYPES.f90
SHARED.o: SHARED.f90 TYPES.o
radiation_module.o: radiation_module.f90 SHARED.o TYPES.o
soil_heat_module.o: soil_heat_module.f90 SHARED.o TYPES.o
energy_balance_module.o: energy_balance_module.f90 SHARED.o TYPES.o radiation_module.o soil_heat_module.o
CANOPY_MODEL.o: CANOPY_MODEL.f90 SHARED.o TYPES.o assim.o th.o meteo.o solar.o energy_balance_module.o
growth_module.o: growth_module.f90 TYPES.o SHARED.o
mineral_module.o: mineral_module.f90 TYPES.o
soiln_module.o: soiln_module.f90 TYPES.o SHARED.o
phenology_module.o: phenology_module.f90 TYPES.o
main_MMCO2Crop.o: main_MMCO2Crop.f90 SHARED.o TYPES.o csv_file.o CANOPY_MODEL.o energy_balance_module.o soil_heat_module.o growth_module.o mineral_module.o soiln_module.o phenology_module.o
th.o: th.f90 SHARED.o TYPES.o
meteo.o: meteo.f90 SHARED.o TYPES.o
solar.o: solar.f90 SHARED.o TYPES.o
hydro.o: hydro.f90 SHARED.o TYPES.o
assim.o: assim.f90 SHARED.o TYPES.o csv_file.o
csv_file.o: csv_file.f90
count_lines.o: count_lines.f90

# Makefile
$(PROG): $(OBJS)
	@echo "------------------------------------"
	@echo " Creating the executable file:      "
	@echo "------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINK)

%.o: %.f90
	@echo "------------------------------------"
	@echo " Compiling file: $<"
	@echo "------------------------------------"
	$(FC) -c $(switch) $(FFLAGS) $<

# Cleaning everything
clean:
	@echo "----------------------------------------"
	@echo " Cleaning garbage files from directory  "
	@echo "----------------------------------------"
	rm -f *.mod *.o $(PROG)
# End of makefile
