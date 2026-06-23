# Compiler
FC = gfortran
FFLAGS = -O3 -Wall

# Directories
SRC = src
TESTS = tests
BUILD = build

# Executables
ROSSBY = rossby.exe
JET = jet.exe

# Default target
all: $(ROSSBY) $(JET)

# Create build directory
$(BUILD):
	mkdir -p $(BUILD)

# -------------------------
# Rossby–Haurwitz test
# -------------------------
$(ROSSBY): $(BUILD)
	$(FC) $(FFLAGS) \
	$(SRC)/solver_mccormack.f90 \
	$(SRC)/grid_yinyang.f90 \
	$(SRC)/stability_strategy1.f90 \
	$(TESTS)/rossby_haurwitz.f90 \
	-o $(ROSSBY)

# -------------------------
# Jet instability test
# -------------------------
$(JET): $(BUILD)
	$(FC) $(FFLAGS) \
	$(SRC)/solver_mccormack.f90 \
	$(SRC)/grid_yinyang.f90 \
	$(SRC)/stability_strategy1.f90 \
	$(TESTS)/jet_instability.f90 \
	-o $(JET)

# Run shortcuts
run-rossby: $(ROSSBY)
	./$(ROSSBY)

run-jet: $(JET)
	./$(JET)

# Clean build
clean:
	rm -f *.exe *.mod
	rm -rf $(BUILD)

# Phony targets
.PHONY: all clean run-rossby run-jet
Add Makefile for build system
