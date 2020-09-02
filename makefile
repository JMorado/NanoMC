#Current make system
# Paths
BIN=./_bin/
SOURCE=./src/
BUILD=./_build/

# The compiler
FC = gfortran
# flags for debugging or for maximum performance  comment as necessary
FCFLAGS = -g -fbounds-check
FCFLAGS = -O2
# flags forall (e.g. look for system .mod files  required in gfortran)
FCFLAGS = -fopenmp

# libraries needed for linking  unused in the examples
#LDFLAGS = 

# List of executables to be built within the package
PROGRAMS = nanomc_nvt.exe nanomc_uvt.exe 
PROGRAMS_PATH=$(addprefix $(BIN), $(PROGRAMS))

# "make" builds all
all: directories $(PROGRAMS_PATH)

# Linking modules
$(BUILD)std_output_module.o: $(BUILD)simulation_module.o
$(BUILD)io_module.o: $(BUILD)simulation_module.o
$(BUILD)pbc_displacement_module.o: $(BUILD)simulation_module.o
$(BUILD)energy_module.o: $(BUILD)pbc_displacement_module.o $(BUILD)simulation_module.o
$(BUILD)monte_carlo_module.o: $(BUILD)simulation_module.o  $(BUILD)energy_module.o  $(BUILD)pbc_displacement_module.o 
$(BUILD)nanomc_nvt.o: $(BUILD)std_output_module.o $(BUILD)io_module.o  $(BUILD)simulation_module.o  $(BUILD)energy_module.o  $(BUILD)pbc_displacement_module.o  $(BUILD)time_module.o  $(BUILD)monte_carlo_module.o
$(BUILD)nanomc_uvt.o: $(BUILD)std_output_module.o $(BUILD)io_module.o  $(BUILD)simulation_module.o  $(BUILD)energy_module.o  $(BUILD)pbc_displacement_module.o  $(BUILD)time_module.o  $(BUILD)monte_carlo_module.o

# Linking modules to executables
$(BIN)nanomc_nvt.exe: $(BUILD)std_output_module.o $(BUILD)io_module.o  $(BUILD)simulation_module.o  $(BUILD)energy_module.o  $(BUILD)pbc_displacement_module.o  $(BUILD)time_module.o  $(BUILD)monte_carlo_module.o
$(BIN)nanomc_uvt.exe: $(BUILD)std_output_module.o $(BUILD)io_module.o  $(BUILD)simulation_module.o  $(BUILD)energy_module.o  $(BUILD)pbc_displacement_module.o  $(BUILD)time_module.o  $(BUILD)monte_carlo_module.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
$(BIN)%.exe: $(BUILD)%.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
$(BUILD)%.o: $(SOURCE)%.f90
	$(FC) $(FCFLAGS) -c $< -o $@ -J $(BUILD)

$(BUILD)%.o: $(SOURCE)%.F90
	$(FC) $(FCFLAGS) -c $< -o $@ -J $(BUILD)

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f $(BUILD)*.o $(BUILD)*.mod $(BUILD)*.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS_PATH)

directories:
	mkdir -p $(BUILD) $(BIN)
