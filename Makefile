# Compiler
FC = gfortran
FFLAGS = -O3 -Wall

# Directories
SRC = src
TESTS = tests

# Output
ADV = advection.exe

# Default build
all: $(ADV)

# -------------------------
# Advection solver
# -------------------------
$(ADV):
	$(FC) $(FFLAGS) \
	$(SRC)/advection_solver.f90 \
	$(SRC)/grid_yinyang.f90 \
	$(TESTS)/advection_test.f90 \
	-o $(ADV)

# Run
run: $(ADV)
	./$(ADV)

# Clean
clean:
	rm -f *.exe *.mod

.PHONY: all run clean
Add Makefile for build system
