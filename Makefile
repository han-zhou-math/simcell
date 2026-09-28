# -*- mode: makefile -*-

# This sample (GNU) Makefile can be used to compile PETSc applications with a single
# source file and can be easily modified to compile multi-file applications.
# It relies on pkg_config tool, and PETSC_DIR and PETSC_ARCH variables.
# Copy this file to your source directory as "Makefile" and modify as needed.
#
# For example - a single source file can be compiled with:
#
#  $ cd src/snes/examples/tutorials/
#  $ make -f $PETSC_DIR/share/petsc/Makefile.user ex17
#
# The following variable must either be a path to PETSc.pc or just "PETSc" if PETSc.pc
# has been installed to a system location or can be found in PKG_CONFIG_PATH.
#PETSC_DIR=/Users/yaol/local/petsc-3.12.5/
#PETSC_ARCH=code-dbg

# =============================================================================
# Compilation instructions (per machine)
# =============================================================================
# PETSC_DIR and PETSC_ARCH must already be set in the environment before `make`
# runs (day-to-day these are exported from the shell profile, not typed by hand).
# The Fortran compiler used here must be the *same one PETSc itself was built
# with* (--with-cc/--with-fc at PETSc configure time) — this code uses modern
# `use petsc`/`use petscsys`/... Fortran modules, and gfortran's .mod file format
# is not compatible across major compiler versions, so a mismatched compiler
# either fails to compile or reads PETSc's .mod files incorrectly.
#
#   macOS (Lingxing's laptop, day-to-day):
#     export PETSC_DIR=/Users/yaol/local/petsc-3.12.5
#     export PETSC_ARCH=osx-gfc-opt
#     make
#
#   OSC Ascend cluster (account PBS0318/uak0394):
#     export PETSC_DIR=/users/PBS0318/uak0394/ESS/local/petsc
#     export PETSC_ARCH=x86_64-gfc-opt
#     module load gcc/12.3.0 openmpi/5.0.2   # matches how PETSc was configured
#                                             # (--with-cc=mpicc --with-fc=mpif90);
#                                             # confirm with `ldd` on libpetsc.so
#                                             # if PETSc is ever rebuilt/updated.
#     make
#
# `make` (= `make all`) runs `clean` then builds `imp` from $(OBJN) below.
# Run with `./imp` (no CLI args, no mpirun needed even though PETSc is built
# with real MPI — reads runtime params from input.par in the cwd).
# =============================================================================

PETSc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/PETSc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(PETSc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
# PETSc 3.19 exposes convergence reasons and vector-array access through the
# legacy integer/F90 interface.  Newer installations use the typed interface.
PETSC_319_CPPFLAGS := $(shell pkg-config --max-version=3.19.99 $(PACKAGES) >/dev/null 2>&1 && echo -DSIMCELL_PETSC_LEGACY_ENUM)
CPPFLAGS += $(PETSC_319_CPPFLAGS)
CPPFLAGS += $(SIMCELL_EXTRA_CPPFLAGS)
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))

print:
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo FC=$(FC)
	@echo CFLAGS=$(CFLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo FFLAGS=$(FFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)

FILES=parameters.f90  myfft.f90 linsys.f90 geometry.f90 chemical_mod.f90 mycontext.f90\
	mycontext_interface.f90 actin_free.f90 network_mod.f90 energy.f90 solver_mod.f90  mainimp.f90 Makefile
#OBJN=parameters.o  myfft.o IBforce.o geometry.o solver_mod.o IBmod.o \
#	linsys.o chemical_mod.o network_mod.o actin_free.o  energy.o mainimp.o
#OBJN=parameters.o  myfft.o mycontext.o mycontext_interface.o geometry.o IBforce.o IBmod.o linsys.o chemical_mod.o fmain.o
# mycontext.f90 / mycontext_interface.f90 are SNES application-context scaffolding for a
# nonlinear solve that was never wired up (FormFunction is an empty stub). Nothing in the
# build uses them any more, so they are kept in the tree for reference but left out of OBJN.
OBJN=parameters.o interface_side_mod.o osmotic_feedback_mod.o fsi_transport_snapshot_mod.o \
	fsi_dualchem_velocity_bridge_mod.o chemical_interface_flux_mod.o \
	chemical_pump_profile_mod.o actin_pnas_profile_mod.o \
	actin_model_mod.o scalar_operator_mod.o generic_gmres_mod.o \
	dualchem/grid_types.o dualchem/small_solver_mod.o dualchem/geometry_mod.o \
	dualchem/advec_diff_solver_mod.o dualchem/linear_solver_mod.o \
	actin_extension_mod.o fsi_actin_sequential_harness_mod.o \
	fsi_actin_block_harness_mod.o fsi_actin_oneway_harness_mod.o \
	fsi_actin_feedback_mod.o \
	fsi_dualchem_coupling_harness_mod.o \
	myfft.o geometry.o IBforce.o IBmod.o linsys.o brinkman_solver_mod.o \
	fsi_actin_force_assembly_mod.o \
	fsi_adhesion_mod.o \
	fsi_external_load_mod.o \
	fsisolve.o fmain.o \
	dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o

# Stage 06 vendors the audited Stage-04 chemical sources under one directory.
# These lists are declared at the import checkpoint for provenance only.  The
# modules are added to build targets after their analytic-velocity dependency
# has been replaced by the tested FSI velocity bridge.
DUALCHEM_DIR := dualchem
DUALCHEM_F90_SRC := $(addprefix $(DUALCHEM_DIR)/,grid_types.f90 \
		small_solver_mod.f90 geometry_mod.f90 advec_diff_solver_mod.f90 \
		linear_solver_mod.f90)
DUALCHEM_CPP_SRC := $(DUALCHEM_DIR)/AdvecDiff2d-semiperiodic/AdvecDiffMG.cpp

FFTWINC=
FFTWDIR=
#FFTWINC=-I${HOME}/local/fftw/gfc/include/
#FFTWDIR=-L${HOME}/local/fftw/gfc/lib/
FFTWINC=-I$(shell pkg-config --variable=includedir fftw3)
FFTWDIR=$(shell pkg-config --libs-only-L fftw3)

#UNAME:=$(shell uname -s)
UNAME:=$(shell uname -s)

ifeq ($(UNAME),Linux)
	ifeq ($(PETSC_ARCH),x86_64-gfc-linux)
		FFLAGS+= -g -O2 -frecursive -fcheck=mem,pointer  -fdefault-real-8\
			 -fdefault-double-8 \
			 -freal-4-real-8 -Wno-unused-variable -Wno-unused-function
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lfftw3 -lmkl_lapack95_lp64 -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
	else ifeq ($(PETSC_ARCH),x86_64-gfc-opt)
		# same as x86_64-gfc-linux but without debug-only checks, and with warnings
		# surfaced instead of silenced, mirroring the osx-gfc-opt convention.
		FFLAGS+= -Wall -Wextra
		# NOTE: plain -llapack/-lblas resolve (via FlexiBLAS on this system) to the
		# ILP64 (64-bit integer) backend, which corrupts args when called from this
		# code's plain 32-bit INTEGER (matches PETSc's own 32-bit PetscInt). Link
		# PETSc's own bundled LP64 fblaslapack instead (already on the -L search
		# path via pkg-config's LDFLAGS, standard unprefixed symbol names).
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lflapack -lfblas -lfftw3
	else ifeq ($(PETSC_ARCH),x86_64-intel-linux)
		FFLAGS+=-O2 -r8 -assume protect_parens  -fp-model strict
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lfftw3 -lmkl_lapack95_lp64 -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
	else
		# Generic pkg-config PETSc installation, including local validation builds.
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lf2clapack -lf2cblas -lfftw3
	endif
endif
ifeq ($(UNAME),Darwin)
#print:
#	@echo "Compiling in OSX..."

	ifeq ($(PETSC_ARCH),osx-gfc)
		FFLAGS+= -g -O2 -frecursive -fcheck=mem,pointer,bounds  -fdefault-real-8 -fdefault-double-8 \
		-freal-4-real-8 -Wno-unused-variable -Wno-unused-function
		FFTWINC=-I${HOME}/local/fftw/gfc/include/
		FFTWDIR=-L${HOME}/local/fftw/gfc/lib/
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lfftw3 -lmkl_lapack95_lp64 -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
	else ifeq ($(PETSC_ARCH),osx-gfc-opt)
		# same as osx-gfc but without debug-only checks, and with warnings surfaced
		# instead of silenced, so gfortran's -Wall/-Wextra actually catch issues.
		FFLAGS+= -Wall -Wextra
		FFTWINC=-I${HOME}/local/fftw/gfc/include/
		FFTWDIR=-L${HOME}/local/fftw/gfc/lib/
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lf2clapack -lf2cblas -lfftw3
	else ifeq ($(PETSC_ARCH),osx_intel)
		FFLAGS+=-O2 -r8 -assume protect_parens  -fp-model strict
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lf2clapack -lf2cblas -lfftw3
#LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lfftw3 -lmkl_lapack95_lp64 -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
		LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm -lmkl_intel_lp64 ##-lf2clapack -lf2cblas -lfftw3
	endif
	#FFLAGS+=-O2 -r8 -assume protect_parens -check bound  -fp-model strict
endif
#FFLAGS+= -g -O2 -frecursive -fcheck=bounds,mem,pointer  -fdefault-real-8 -fdefault-double-8 \
#	-freal-4-real-8 -Wno-unused-variable -Wno-unused-function
#FFLAGS+=-O2 -r8 -assume protect_parens -check bound  -fp-model strict

.PHONY: all clean test_shared_parameters test_geometry_origin \
	check_pnas_pump_profile check_pnas_actin_profile check_interface_side \
	check_chemical_diffusion_scaling \
	test_stokes_fourier_phase test_stokes_boundary_extension \
	check_fsi_dualchem_velocity_bridge check_dualchem_velocity_routing \
	check_fsi_dualchem_oneway check_fsi_transport_snapshot \
	check_stage09_shared_snapshot check_stage09_acceptance_tests \
	check_actin_model check_scalar_operator_contract check_generic_gmres \
	check_actin_local_correction \
	check_actin_scalar_extension check_fsi_actin_sequential \
	check_actin_reaction_transfer check_actin_static_exchange \
	check_actin_moving_history \
	check_fsi_actin_block \
	check_fsi_actin_oneway_harness \
	check_stage12_shared_snapshot \
	check_stage12_analysis_tests \
	check_stage12_refinement_analysis_tests \
	check_brinkman_operator \
	check_brinkman_lagged_comparison \
	check_stage13_routing_tests \
	check_stage13_zero_drag_analysis_tests \
	check_fsi_brinkman_routing \
	check_fsi_adhesion check_fsi_external_load \
	check_actin_fsi_feedback_builder \
	check_actin_fsi_force_assembly \
	check_stage14_feedback_routing_tests \
	check_stage14_feedback_transaction \
	check_stage14_analysis_tests \
	dualchem_modules

all: clean imp

imp: $(OBJN)
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS) \
		$(STAGE06_CXX_RUNTIME)

#$(COMPILE.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) -lfftw3 $(LDLIBS)

clean:
	@rm -f $(OBJS) imp pc test_shared_parameters test_geometry_origin test_valatibpt_fixed_stencil \
		test_interface_side test_pnas_pump_profile test_pnas_actin_profile \
		test_chemical_diffusion_scaling \
		test_stokes_fourier_phase test_stokes_boundary_extension \
	test_fsi_dualchem_velocity_bridge test_dualchem_velocity_routing \
	test_fsi_dualchem_oneway test_physical_concentration_jump \
	test_fsi_transport_snapshot test_stage09_shared_snapshot \
	test_osmotic_slip_sign test_fsi_dualchem_twoway \
		test_actin_model \
		test_scalar_operator_contract \
		test_generic_gmres \
		test_actin_local_correction \
		test_actin_scalar_extension \
		test_fsi_actin_sequential \
		test_actin_reaction_transfer \
		test_actin_static_exchange \
		test_actin_moving_history \
		test_fsi_actin_block \
		test_fsi_actin_oneway_harness \
		test_stage12_shared_snapshot \
		test_brinkman_operator \
		test_brinkman_lagged_comparison \
		test_fsi_brinkman_routing \
		test_fsi_adhesion \
		test_fsi_external_load \
		test_actin_fsi_feedback_builder \
		test_actin_fsi_force_assembly \
		test_stage14_feedback_transaction \
		*.dat dump.m *.mod *.o dualchem/*.o \
		dualchem/AdvecDiff2d-semiperiodic/*.o *.s log
	@rm -rf tests/__pycache__

myfft.o: myfft.f90
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp $(FFTWINC)

fsi_adhesion_mod.o: fsi_adhesion_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

fsi_external_load_mod.o: fsi_external_load_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

interface_side_mod.o: interface_side_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_interface_side.o: test_interface_side.f90 parameters.o interface_side_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_interface_side: parameters.o interface_side_mod.o test_interface_side.o
	$(FC) $(FFLAGS) -o $@ $^

check_interface_side: test_interface_side
	./test_interface_side

test_actin_fsi_feedback_builder test_actin_fsi_force_assembly \
	test_actin_local_correction test_actin_moving_history \
	test_actin_reaction_transfer test_actin_scalar_extension \
	test_actin_static_exchange test_brinkman_lagged_comparison \
	test_brinkman_operator test_dualchem_velocity_routing \
	test_fsi_actin_block test_fsi_actin_oneway_harness \
	test_fsi_actin_sequential test_fsi_brinkman_routing \
	test_fsi_dualchem_oneway test_fsi_dualchem_twoway \
	test_geometry_origin test_physical_concentration_jump \
	test_stage09_shared_snapshot test_stage12_shared_snapshot \
	test_stage14_feedback_transaction: interface_side_mod.o

test_fsi_adhesion.o: test_fsi_adhesion.f90 parameters.o fsi_adhesion_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_adhesion: parameters.o fsi_adhesion_mod.o test_fsi_adhesion.o
	$(FC) $(FFLAGS) -o $@ $^

check_fsi_adhesion: test_fsi_adhesion
	./test_fsi_adhesion

test_fsi_external_load.o: test_fsi_external_load.f90 parameters.o \
		fsi_external_load_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_external_load: parameters.o fsi_external_load_mod.o \
		test_fsi_external_load.o
	$(FC) $(FFLAGS) -o $@ $^

check_fsi_external_load: test_fsi_external_load
	./test_fsi_external_load

check_stage09_acceptance_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
		tests/test_check_stage09_acceptance.py -v

# Stage-06 checkpoint: verify the MAC-layout map and marker interpolation before
# any dualchem production routine is allowed to consume the FSI velocity.
fsi_dualchem_velocity_bridge_mod.o: fsi_dualchem_velocity_bridge_mod.f90 \
		parameters.o fsi_transport_snapshot_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

chemical_interface_flux_mod.o: chemical_interface_flux_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_chemical_diffusion_scaling.o: test_chemical_diffusion_scaling.f90 \
		parameters.o chemical_interface_flux_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_chemical_diffusion_scaling: parameters.o chemical_interface_flux_mod.o \
		test_chemical_diffusion_scaling.o
	$(FC) $(FFLAGS) -o $@ $^

check_chemical_diffusion_scaling: test_chemical_diffusion_scaling
	./test_chemical_diffusion_scaling

chemical_pump_profile_mod.o: chemical_pump_profile_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_pnas_pump_profile.o: test_pnas_pump_profile.f90 parameters.o \
		chemical_pump_profile_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_pnas_pump_profile: parameters.o chemical_pump_profile_mod.o \
		test_pnas_pump_profile.o
	$(FC) $(FFLAGS) -o $@ $^

check_pnas_pump_profile: test_pnas_pump_profile
	./test_pnas_pump_profile

actin_pnas_profile_mod.o: actin_pnas_profile_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_pnas_actin_profile.o: test_pnas_actin_profile.f90 parameters.o \
		actin_pnas_profile_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_pnas_actin_profile: parameters.o actin_pnas_profile_mod.o \
		test_pnas_actin_profile.o
	$(FC) $(FFLAGS) -o $@ $^

check_pnas_actin_profile: test_pnas_actin_profile
	./test_pnas_actin_profile

fsi_transport_snapshot_mod.o: fsi_transport_snapshot_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_transport_snapshot.o: test_fsi_transport_snapshot.f90 parameters.o \
		fsi_transport_snapshot_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_transport_snapshot: parameters.o fsi_transport_snapshot_mod.o \
		test_fsi_transport_snapshot.o
	$(FC) $(FFLAGS) -o $@ $^

check_fsi_transport_snapshot: test_fsi_transport_snapshot
	./test_fsi_transport_snapshot

test_stage09_shared_snapshot.o: test_stage09_shared_snapshot.f90 parameters.o \
		fsi_transport_snapshot_mod.o fsi_dualchem_velocity_bridge_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_stage09_shared_snapshot: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o test_stage09_shared_snapshot.o
	$(FC) $(FFLAGS) -o $@ $^

check_stage09_shared_snapshot: test_stage09_shared_snapshot
	./test_stage09_shared_snapshot

# Stage 10 checkpoint A: reduced actin coefficients and the exact physical/
# normalized moving-interface flux signs, independent of every PDE backend.
actin_model_mod.o: actin_model_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_model.o: test_actin_model.f90 parameters.o actin_model_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_model: parameters.o actin_model_mod.o test_actin_model.o
	$(FC) $(FFLAGS) -o $@ $^

check_actin_model: test_actin_model
	./test_actin_model

scalar_operator_mod.o: scalar_operator_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_scalar_operator_contract.o: test_scalar_operator_contract.f90 \
		parameters.o scalar_operator_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_scalar_operator_contract: parameters.o scalar_operator_mod.o \
		test_scalar_operator_contract.o
	$(FC) $(FFLAGS) -o $@ $^

check_scalar_operator_contract: test_scalar_operator_contract
	./test_scalar_operator_contract

generic_gmres_mod.o: generic_gmres_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_generic_gmres.o: test_generic_gmres.f90 parameters.o generic_gmres_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_generic_gmres: parameters.o generic_gmres_mod.o test_generic_gmres.o
	$(FC) $(FFLAGS) -o $@ $^

check_generic_gmres: test_generic_gmres
	./test_generic_gmres

.PHONY: check_valatibpt_fixed_stencil
check_valatibpt_fixed_stencil: test_valatibpt_fixed_stencil
	./test_valatibpt_fixed_stencil

test_valatibpt_fixed_stencil.o: test_valatibpt_fixed_stencil.f90 \
		parameters.o dualchem/grid_types.o dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_valatibpt_fixed_stencil: parameters.o interface_side_mod.o \
		fsi_transport_snapshot_mod.o fsi_dualchem_velocity_bridge_mod.o \
		chemical_interface_flux_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		test_valatibpt_fixed_stencil.o
	$(FC) $(FFLAGS) -o $@ $^

test_actin_local_correction.o: test_actin_local_correction.f90 \
		parameters.o dualchem/grid_types.o dualchem/small_solver_mod.o \
		dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_local_correction: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		test_actin_local_correction.o
	$(FC) $(FFLAGS) -o $@ $^

check_actin_local_correction: test_actin_local_correction
	./test_actin_local_correction

actin_extension_mod.o: actin_extension_mod.f90 parameters.o \
		scalar_operator_mod.o dualchem/grid_types.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_scalar_extension.o: test_actin_scalar_extension.f90 parameters.o \
		scalar_operator_mod.o actin_extension_mod.o dualchem/grid_types.o \
		dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_scalar_extension: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o scalar_operator_mod.o \
		dualchem/grid_types.o dualchem/small_solver_mod.o \
		dualchem/geometry_mod.o dualchem/advec_diff_solver_mod.o \
		actin_extension_mod.o dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_actin_scalar_extension.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_actin_scalar_extension: test_actin_scalar_extension
	./test_actin_scalar_extension

# Stage 10 checkpoint E: network/F-actin is solved and explicitly finalized
# before free/G-actin sees gamma times that staged network field.  The focused
# test uses an injected affine extension to make ordering and rollback exact.
fsi_actin_sequential_harness_mod.o: fsi_actin_sequential_harness_mod.f90 \
		parameters.o fsi_transport_snapshot_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o scalar_operator_mod.o generic_gmres_mod.o \
		actin_extension_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_sequential.o: test_fsi_actin_sequential.f90 parameters.o \
		fsi_transport_snapshot_mod.o actin_model_mod.o scalar_operator_mod.o \
		generic_gmres_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_sequential: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_fsi_actin_sequential.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_fsi_actin_sequential: test_fsi_actin_sequential
	./test_fsi_actin_sequential

# Stage 10 physical validation: the approved network-first/current-network
# reaction split conserves theta_n+theta_c exactly in its discrete algebra.
test_actin_reaction_transfer.o: test_actin_reaction_transfer.f90 parameters.o \
		fsi_transport_snapshot_mod.o actin_model_mod.o scalar_operator_mod.o \
		generic_gmres_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		dualchem/grid_types.o \
		dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_reaction_transfer: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_actin_reaction_transfer.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_actin_reaction_transfer: test_actin_reaction_transfer
	./test_actin_reaction_transfer

# Stage 10 physical validation: a static circular interface exercises the two
# Robin laws through the real Fortran-to-C++ semi-periodic scalar backend.
test_actin_static_exchange.o: test_actin_static_exchange.f90 parameters.o \
		fsi_transport_snapshot_mod.o actin_model_mod.o scalar_operator_mod.o \
		generic_gmres_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		dualchem/grid_types.o \
		dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_static_exchange: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_actin_static_exchange.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_actin_static_exchange: test_actin_static_exchange
	./test_actin_static_exchange

# Stage 10 moving-history validation starts with the two swept-cell signs and
# the frozen old nearest-marker/old-polynomial-center contract.
test_actin_moving_history.o: test_actin_moving_history.f90 parameters.o \
		dualchem/grid_types.o dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_moving_history: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		test_actin_moving_history.o
	$(FC) $(FFLAGS) -o $@ $^ -lf2clapack -lf2cblas

check_actin_moving_history: test_actin_moving_history
	./test_actin_moving_history

# Stage 11 adds only the packed boundary-density harness.  It reuses the
# Stage 10 species evaluator and unchanged semi-periodic Cartesian backend.
fsi_actin_block_harness_mod.o: fsi_actin_block_harness_mod.f90 parameters.o \
		actin_model_mod.o scalar_operator_mod.o generic_gmres_mod.o \
		actin_extension_mod.o fsi_actin_sequential_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_block.o: test_fsi_actin_block.f90 parameters.o \
		fsi_actin_block_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_block: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_fsi_actin_block.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_fsi_actin_block: test_fsi_actin_block
	./test_fsi_actin_block

# Stage 12 persistent one-way manager.  The manager owns only accepted actin
# state/geometry and consumes the snapshot already accepted by dualchem.
fsi_actin_oneway_harness_mod.o: fsi_actin_oneway_harness_mod.f90 parameters.o \
		fsi_transport_snapshot_mod.o actin_model_mod.o dualchem/grid_types.o \
		dualchem/geometry_mod.o fsi_actin_sequential_harness_mod.o \
		fsi_actin_block_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_oneway_harness.o: test_fsi_actin_oneway_harness.f90 \
		parameters.o fsi_actin_oneway_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_actin_oneway_harness: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		fsi_actin_oneway_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_fsi_actin_oneway_harness.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_fsi_actin_oneway_harness: test_fsi_actin_oneway_harness
	./test_fsi_actin_oneway_harness

# Stage 12 integration gate: the chemical transaction remains the sole
# snapshot publisher and the persistent actin manager consumes that exact id.
test_stage12_shared_snapshot.o: test_stage12_shared_snapshot.f90 parameters.o \
		fsi_dualchem_coupling_harness_mod.o fsi_actin_oneway_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_stage12_shared_snapshot: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o dualchem/linear_solver_mod.o \
		chemical_interface_flux_mod.o \
		fsi_dualchem_coupling_harness_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		fsi_actin_oneway_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_stage12_shared_snapshot.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_stage12_shared_snapshot: test_stage12_shared_snapshot
	./test_stage12_shared_snapshot

check_stage12_analysis_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_check_stage12_oneway.py' -v

check_stage12_refinement_analysis_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_analyze_stage12_refinement.py' -v

# Stage 13 solves the variable-coefficient Brinkman problem in the
# divergence-free velocity space while reusing the accepted rectangular
# Fourier/banded Stokes inverse for every projected operator application.
brinkman_solver_mod.o: brinkman_solver_mod.f90 parameters.o scalar_operator_mod.o \
		myfft.o IBforce.o IBmod.o linsys.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_brinkman_operator.o: test_brinkman_operator.f90 parameters.o \
		brinkman_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_brinkman_operator: parameters.o scalar_operator_mod.o myfft.o IBforce.o \
		IBmod.o linsys.o brinkman_solver_mod.o test_brinkman_operator.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

check_brinkman_operator: test_brinkman_operator
	./test_brinkman_operator

test_brinkman_lagged_comparison.o: test_brinkman_lagged_comparison.f90 \
		parameters.o brinkman_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_brinkman_lagged_comparison: parameters.o scalar_operator_mod.o myfft.o \
		IBforce.o IBmod.o linsys.o brinkman_solver_mod.o \
		test_brinkman_lagged_comparison.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

check_brinkman_lagged_comparison: test_brinkman_lagged_comparison
	./test_brinkman_lagged_comparison

check_stage13_routing_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_stage13_brinkman_routing.py' -v

check_stage13_zero_drag_analysis_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_compare_stage13_zero_drag.py' -v

test_fsi_brinkman_routing.o: test_fsi_brinkman_routing.f90 parameters.o \
		brinkman_solver_mod.o fsi_actin_feedback_mod.o fsisolve.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_brinkman_routing: parameters.o scalar_operator_mod.o myfft.o \
		IBforce.o IBmod.o linsys.o osmotic_feedback_mod.o \
		brinkman_solver_mod.o actin_model_mod.o fsi_actin_feedback_mod.o \
		fsi_actin_force_assembly_mod.o fsisolve.o test_fsi_brinkman_routing.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

check_fsi_brinkman_routing: test_fsi_brinkman_routing
	./test_fsi_brinkman_routing

# Stage 14 checkpoint A: export the accepted network state and construct the
# frozen MAC-face Brinkman coefficient, bulk stress-gradient force, and marker
# stress.  No FSI solve or scalar C++ backend is changed at this checkpoint.
test_actin_fsi_feedback_builder.o: test_actin_fsi_feedback_builder.f90 \
		parameters.o actin_model_mod.o fsi_actin_oneway_harness_mod.o \
		fsi_actin_feedback_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

fsi_actin_feedback_mod.o: fsi_actin_feedback_mod.f90 parameters.o \
		scalar_operator_mod.o actin_model_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_fsi_feedback_builder: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o actin_model_mod.o \
		actin_pnas_profile_mod.o \
		scalar_operator_mod.o generic_gmres_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o actin_extension_mod.o \
		fsi_actin_sequential_harness_mod.o fsi_actin_block_harness_mod.o \
		fsi_actin_oneway_harness_mod.o fsi_actin_feedback_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_actin_fsi_feedback_builder.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_actin_fsi_feedback_builder: test_actin_fsi_feedback_builder
	./test_actin_fsi_feedback_builder

test_actin_fsi_force_assembly.o: test_actin_fsi_force_assembly.f90 \
		parameters.o myfft.o IBforce.o IBmod.o actin_model_mod.o \
		fsi_actin_feedback_mod.o fsi_actin_force_assembly_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

fsi_actin_force_assembly_mod.o: fsi_actin_force_assembly_mod.f90 \
		parameters.o myfft.o IBforce.o IBmod.o scalar_operator_mod.o \
		fsi_actin_feedback_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_actin_fsi_force_assembly: parameters.o myfft.o IBforce.o IBmod.o \
		scalar_operator_mod.o actin_model_mod.o fsi_actin_feedback_mod.o \
		fsi_actin_force_assembly_mod.o \
		test_actin_fsi_force_assembly.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

check_actin_fsi_force_assembly: test_actin_fsi_force_assembly
	./test_actin_fsi_force_assembly

check_stage14_feedback_routing_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_stage14_feedback_routing.py' -v

# Stage 14 checkpoint D: run identical one-step implicit FSI transactions with
# disabled and active accepted feedback, then repeat the active case.  This
# exercises PETSc ordinary/MatShell routing, source-id retention, rollback at
# the pre-FSI build boundary, and deterministic publication.
test_stage14_feedback_transaction.o: test_stage14_feedback_transaction.f90 \
		parameters.o actin_model_mod.o fsi_actin_feedback_mod.o fsisolve.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_stage14_feedback_transaction: parameters.o scalar_operator_mod.o myfft.o \
		IBforce.o IBmod.o linsys.o osmotic_feedback_mod.o \
		brinkman_solver_mod.o actin_model_mod.o fsi_actin_feedback_mod.o \
		fsi_actin_force_assembly_mod.o fsisolve.o \
		test_stage14_feedback_transaction.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

check_stage14_feedback_transaction: test_stage14_feedback_transaction
	./test_stage14_feedback_transaction
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_stage14_feedback_repeat.py' -v

check_stage14_analysis_tests:
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_check_stage14_feedback.py' -v
	PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests \
		-p 'test_stage14_study_analysis.py' -v

osmotic_feedback_mod.o: osmotic_feedback_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

fsisolve.o: fsisolve.f90 osmotic_feedback_mod.o brinkman_solver_mod.o \
		fsi_actin_feedback_mod.o fsi_actin_force_assembly_mod.o \
		fsi_adhesion_mod.o fsi_external_load_mod.o

test_fsi_dualchem_velocity_bridge.o: test_fsi_dualchem_velocity_bridge.f90 \
		parameters.o fsi_dualchem_velocity_bridge_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_dualchem_velocity_bridge: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o \
		test_fsi_dualchem_velocity_bridge.o
	$(FC) $(FFLAGS) -o $@ $^

check_fsi_dualchem_velocity_bridge: test_fsi_dualchem_velocity_bridge
	./test_fsi_dualchem_velocity_bridge

dualchem/grid_types.o: dualchem/grid_types.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

dualchem/small_solver_mod.o: dualchem/small_solver_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

dualchem/geometry_mod.o: dualchem/geometry_mod.f90 parameters.o \
		chemical_interface_flux_mod.o interface_side_mod.o \
		fsi_dualchem_velocity_bridge_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

geometry.o: geometry.f90 parameters.o interface_side_mod.o myfft.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

dualchem/advec_diff_solver_mod.o: dualchem/advec_diff_solver_mod.f90 parameters.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

dualchem/linear_solver_mod.o: dualchem/linear_solver_mod.f90 parameters.o \
		chemical_interface_flux_mod.o chemical_pump_profile_mod.o \
		fsi_dualchem_velocity_bridge_mod.o \
		dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o: \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.cpp
	$(CXX) $(CXXFLAGS) -Idualchem/AdvecDiff2d-semiperiodic/include -c -o $@ $<

dualchem_modules: dualchem/grid_types.o dualchem/small_solver_mod.o \
		dualchem/geometry_mod.o dualchem/advec_diff_solver_mod.o \
		dualchem/linear_solver_mod.o chemical_interface_flux_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o

test_dualchem_velocity_routing.o: test_dualchem_velocity_routing.f90 parameters.o \
		fsi_dualchem_velocity_bridge_mod.o dualchem/grid_types.o dualchem/geometry_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_dualchem_velocity_routing: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o \
		dualchem/grid_types.o dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		test_dualchem_velocity_routing.o
	$(FC) $(FFLAGS) -o $@ $^

check_dualchem_velocity_routing: test_dualchem_velocity_routing
	./test_dualchem_velocity_routing

fsi_dualchem_coupling_harness_mod.o: fsi_dualchem_coupling_harness_mod.f90 \
		parameters.o fsi_transport_snapshot_mod.o \
		chemical_interface_flux_mod.o fsi_dualchem_velocity_bridge_mod.o \
		dualchem/grid_types.o \
		dualchem/geometry_mod.o dualchem/linear_solver_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

fmain.o: fmain.f90 fsi_dualchem_coupling_harness_mod.o \
		fsi_actin_oneway_harness_mod.o actin_model_mod.o \
		fsi_actin_feedback_mod.o osmotic_feedback_mod.o \
		brinkman_solver_mod.o fsisolve.o

test_fsi_dualchem_oneway.o: test_fsi_dualchem_oneway.f90 parameters.o \
		fsi_dualchem_coupling_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

ifeq ($(UNAME),Darwin)
STAGE06_CXX_RUNTIME=-lc++
else
STAGE06_CXX_RUNTIME=-lstdc++
endif

test_fsi_dualchem_oneway: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o \
		dualchem/grid_types.o dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o dualchem/linear_solver_mod.o \
		chemical_interface_flux_mod.o \
		fsi_dualchem_coupling_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_fsi_dualchem_oneway.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_fsi_dualchem_oneway: test_fsi_dualchem_oneway
	./test_fsi_dualchem_oneway

test_physical_concentration_jump.o: test_physical_concentration_jump.f90 \
		parameters.o fsi_dualchem_coupling_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_physical_concentration_jump: parameters.o fsi_transport_snapshot_mod.o \
		fsi_dualchem_velocity_bridge_mod.o dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o dualchem/linear_solver_mod.o \
		chemical_interface_flux_mod.o \
		fsi_dualchem_coupling_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_physical_concentration_jump.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_physical_concentration_jump: test_physical_concentration_jump
	./test_physical_concentration_jump

test_osmotic_slip_sign.o: test_osmotic_slip_sign.f90 parameters.o \
		osmotic_feedback_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_osmotic_slip_sign: parameters.o osmotic_feedback_mod.o \
		test_osmotic_slip_sign.o
	$(FC) $(FFLAGS) -o $@ $^

check_osmotic_slip_sign: test_osmotic_slip_sign
	./test_osmotic_slip_sign

test_fsi_dualchem_twoway.o: test_fsi_dualchem_twoway.f90 parameters.o \
		osmotic_feedback_mod.o fsi_dualchem_coupling_harness_mod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_fsi_dualchem_twoway: parameters.o osmotic_feedback_mod.o \
		fsi_transport_snapshot_mod.o fsi_dualchem_velocity_bridge_mod.o \
		dualchem/grid_types.o \
		dualchem/small_solver_mod.o dualchem/geometry_mod.o \
		dualchem/advec_diff_solver_mod.o dualchem/linear_solver_mod.o \
		chemical_interface_flux_mod.o \
		fsi_dualchem_coupling_harness_mod.o \
		dualchem/AdvecDiff2d-semiperiodic/AdvecDiffMG.o \
		test_fsi_dualchem_twoway.o
	$(FC) $(FFLAGS) -o $@ $^ $(STAGE06_CXX_RUNTIME)

check_fsi_dualchem_twoway: test_fsi_dualchem_twoway
	./test_fsi_dualchem_twoway

test_shared_parameters.o: test_shared_parameters.f90 parameters.o IBmod.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_shared_parameters: parameters.o myfft.o IBforce.o IBmod.o test_shared_parameters.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

# This focused executable exercises the actual geometry linked-list setup.
# It deliberately poisons the legacy xamin/yamin aliases before getIBlist:
# grid indices must be derived from the immutable domain origin xmin/ymin.
test_geometry_origin.o: test_geometry_origin.f90 parameters.o geometry.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_geometry_origin: parameters.o myfft.o geometry.o IBforce.o IBmod.o \
		linsys.o osmotic_feedback_mod.o fsisolve.o test_geometry_origin.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

# The Fourier phase is a grid-index phase, 2*pi*k/nx, and must not change
# when the physical periodic length changes from one to two.
test_stokes_fourier_phase.o: test_stokes_fourier_phase.f90 parameters.o linsys.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_stokes_fourier_phase: parameters.o myfft.o IBforce.o IBmod.o linsys.o \
		test_stokes_fourier_phase.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)

test_stokes_boundary_extension.o: test_stokes_boundary_extension.f90 parameters.o linsys.o
	$(COMPILE.F) $(OUTPUT_OPTION) $< -cpp

test_stokes_boundary_extension: parameters.o myfft.o IBforce.o IBmod.o linsys.o \
		test_stokes_boundary_extension.o
	$(LINK.F) $(FFLAGS) -o $@ $^ $(FFTWINC) $(FFTWDIR) $(LDLIBS)



# Many suffixes are covered by implicit rules, but you may need to write custom rules
# such as these if you use suffixes that do not have implicit rules.
# https://www.gnu.org/software/make/manual/html_node/Catalogue-of-Rules.html#Catalogue-of-Rules


% : %.f90
	$(LINK.F) -o $@ $^ $(LDLIBS)
%.o: %.f90
	$(COMPILE.F) $(OUTPUT_OPTION) $<  -cpp
% : %.cxx
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cxx
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

# For a multi-file case, suppose you have the source files a.F90, b.c, and c.cxx
# (with a program statement appearing in a.F90 or main() appearing in the C or
# C++ source).  This can be built by uncommenting the following two lines.
#
# app : a.o b.o c.o
# 	$(LINK.F) -o $@ $^ $(LDLIBS)

# If the file c.cxx needs to link with a C++ standard library -lstdc++ , then
# you'll need to add it explicitly.  It can go in the rule above or be added to
# a target-specific variable by uncommenting the line below.
#
# app : LDLIBS += -lstdc++
