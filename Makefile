# Makefile for Fortran project
# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -Wall -g -O2

# Source files
SOURCES = mod_var.f90 mod_function.f90 mod_fem.f90 main_v2.f90

# Object files (replace .f90 with .o)
OBJECTS = $(SOURCES:.f90=.o)

# Module files (will be created during compilation)
MODULES = mod_var.mod mod_function.mod mod_fem.mod

# Executable name
EXECUTABLE = wave_solver

# Default target
all: $(EXECUTABLE)

# Build executable
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Generic rule for object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Dependencies (modules must be compiled before files that use them)
mod_function.o: mod_var.o
mod_fem.o: mod_var.o mod_function.o
main_v2.o: mod_var.o mod_function.o mod_fem.o

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(MODULES) $(EXECUTABLE)

# Clean and rebuild
rebuild: clean all

# Show help
help:
	@echo "Available targets:"
	@echo "  all      - Build the executable (default)"
	@echo "  clean    - Remove object files, modules, and executable"
	@echo "  rebuild  - Clean and build"
	@echo "  help     - Show this help message"

# Phony targets (not files)
.PHONY: all clean rebuild help