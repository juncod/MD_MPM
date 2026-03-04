# Makefile for LAMMPS-MPM Electric Field Coupling
#
# Dependencies:
#   - LAMMPS (shared library): set LAMMPS_DIR
#   - SuiteSparse (UMFPACK): set SUITESPARSE_DIR
#   - MPI: mpicc
#
# Usage:
#   make all          — build main executable
#   make test         — build and run unit tests
#   make test_mesh    — mesh unit test only
#   make test_basis   — basis function unit test only
#   make test_solver  — solver integration test only
#   make clean        — remove build artifacts

# ─── Paths (adjust for your system) ───────────────────────────────
LAMMPS_DIR      ?= /usr/local
SUITESPARSE_DIR ?= /usr/local

# ─── Compilers ────────────────────────────────────────────────────
CC      = mpicc
CFLAGS  = -Wall -O2 -std=c99
LDFLAGS =

# Test compiler (no MPI needed for mesh/basis tests)
TCC     = cc
TCFLAGS = -Wall -O2 -std=c99

# ─── Include / Library paths ─────────────────────────────────────
INCLUDES = -I./include \
           -I$(LAMMPS_DIR)/include \
           -I$(SUITESPARSE_DIR)/include \
           -I$(SUITESPARSE_DIR)/include/suitesparse

LIBS_UMFPACK = -L$(SUITESPARSE_DIR)/lib \
               -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd \
               -lsuitesparseconfig

LIBS_LAMMPS = -L$(LAMMPS_DIR)/lib -llammps

LIBS = $(LIBS_UMFPACK) $(LIBS_LAMMPS) -lm

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

# ─── Unit Tests ──────────────────────────────────────────────────

# Mesh test (no UMFPACK, no MPI, no LAMMPS)
test_mesh: dirs $(TEST_DIR)/test_mesh.c $(SRC_DIR)/mpm_mesh.c $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include -o $(BUILD_DIR)/test_mesh \
		$(TEST_DIR)/test_mesh.c $(SRC_DIR)/mpm_mesh.c -lm
	@echo "Running mesh tests..."
	$(BUILD_DIR)/test_mesh

# Basis test (no UMFPACK, no MPI, no LAMMPS)
test_basis: dirs $(TEST_DIR)/test_basis.c $(SRC_DIR)/mpm_basis.c $(SRC_DIR)/mpm_mesh.c $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include -o $(BUILD_DIR)/test_basis \
		$(TEST_DIR)/test_basis.c $(SRC_DIR)/mpm_basis.c $(SRC_DIR)/mpm_mesh.c -lm
	@echo "Running basis tests..."
	$(BUILD_DIR)/test_basis

# Solver test (needs UMFPACK, no MPI/LAMMPS)
test_solver: dirs $(TEST_DIR)/test_solver.c $(SRCS_CORE) $(HEADERS)
	$(TCC) $(TCFLAGS) -I./include $(INCLUDES) -o $(BUILD_DIR)/test_solver \
		$(TEST_DIR)/test_solver.c $(SRCS_CORE) $(LIBS_UMFPACK) -lm
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
