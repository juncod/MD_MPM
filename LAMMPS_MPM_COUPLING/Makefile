# Makefile for LAMMPS-MPM Electric Field Coupling
#
# Dependencies:
#   - LAMMPS (shared library): set LAMMPS_DIR
#   - MPI: mpicc
#   - No external solver library needed (built-in CG solver)
#
# Usage:
#   make all          — build main executable
#   make test         — build and run unit tests
#   make test_mesh    — mesh unit test only
#   make test_basis   — basis function unit test only
#   make test_solver  — solver integration test only
#   make clean        — remove build artifacts

# ─── Paths (adjust for your system) ───────────────────────────────
LAMMPS_DIR ?= /usr/local

# ─── Compilers ────────────────────────────────────────────────────
CC      = mpicc
CFLAGS  = -Wall -O2 -std=c99
LDFLAGS =

# Test compiler (no MPI needed for unit tests)
TCC     = cc
TCFLAGS = -Wall -O2 -std=c99

# ─── Include / Library paths ─────────────────────────────────────
INCLUDES = -I./include \
           -I$(LAMMPS_DIR)/include

LIBS_LAMMPS = -L$(LAMMPS_DIR)/lib -llammps

LIBS = $(LIBS_LAMMPS) -lm

# ─── Sources ─────────────────────────────────────────────────────
SRC_DIR   = src
INC_DIR   = include
TEST_DIR  = tests
BUILD_DIR = build

SRCS = $(SRC_DIR)/main.c \
       $(SRC_DIR)/mpm_mesh.c \
       $(SRC_DIR)/mpm_basis.c \
       $(SRC_DIR)/mpm_assembly.c \
       $(SRC_DIR)/mpm_solver.c \
       $(SRC_DIR)/mpm_electric.c \
       $(SRC_DIR)/lammps_callback.c \
       $(SRC_DIR)/atom_to_mp.c \
       $(SRC_DIR)/vtk_output.c

# Core MPM sources (no LAMMPS/MPI dependency)
SRCS_CORE = $(SRC_DIR)/mpm_mesh.c \
            $(SRC_DIR)/mpm_basis.c \
            $(SRC_DIR)/mpm_assembly.c \
            $(SRC_DIR)/mpm_solver.c \
            $(SRC_DIR)/mpm_electric.c \
            $(SRC_DIR)/atom_to_mp.c \
            $(SRC_DIR)/vtk_output.c

HEADERS = $(wildcard $(INC_DIR)/*.h)

# ─── Targets ─────────────────────────────────────────────────────
TARGET = lammps_mpm_coupling

.PHONY: all clean test test_mesh test_basis test_solver dirs

all: dirs $(TARGET)

dirs:
	@mkdir -p $(BUILD_DIR) output

$(TARGET): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BUILD_DIR)/$(TARGET) $(SRCS) $(LDFLAGS) $(LIBS)
	@echo "Built: $(BUILD_DIR)/$(TARGET)"

# ─── Unit Tests (no external dependencies except -lm) ────────────

# Mesh test
test_mesh: dirs $(TEST_DIR)/test_mesh.c $(SRC_DIR)/mpm_mesh.c $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include -o $(BUILD_DIR)/test_mesh \
		$(TEST_DIR)/test_mesh.c $(SRC_DIR)/mpm_mesh.c -lm
	@echo "Running mesh tests..."
	$(BUILD_DIR)/test_mesh

# Basis test
test_basis: dirs $(TEST_DIR)/test_basis.c $(SRC_DIR)/mpm_basis.c $(SRC_DIR)/mpm_mesh.c $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include -o $(BUILD_DIR)/test_basis \
		$(TEST_DIR)/test_basis.c $(SRC_DIR)/mpm_basis.c $(SRC_DIR)/mpm_mesh.c -lm
	@echo "Running basis tests..."
	$(BUILD_DIR)/test_basis

# Solver integration test (no external libs needed — pure C CG solver)
test_solver: dirs $(TEST_DIR)/test_solver.c $(SRCS_CORE) $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include -o $(BUILD_DIR)/test_solver \
		$(TEST_DIR)/test_solver.c $(SRCS_CORE) -lm
	@echo "Running solver tests..."
	$(BUILD_DIR)/test_solver

# Run all tests
test: test_mesh test_basis test_solver
	@echo ""
	@echo "=== All tests passed ==="

# ─── Clean ───────────────────────────────────────────────────────
clean:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned build artifacts"
