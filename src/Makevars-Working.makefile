# --- 3. Custom GNU Make Dependency Rules ---
# These rules force the compilation order for modules that are used by others.
# Note: R's build system typically handles the final linking step, but we must
# ensure modules are built and used in the correct order using these dependencies.

.PHONY: all
all: $(SHLIB)

# --- Module Compilation Rules ---
# These rules explicitly define how the module file (.mod) and object file (.o)
# are created, and on which other modules they depend.

# 0. Base Parameters Module (NEW RULE)
# This module must be built first as all other modules depend on it.
rprintf_mod.o rprintf_mod.mod: rprintf_mod.f90
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o rprintf_mod.o

00tweedie_params.o 00tweedie_params.mod: 00tweedie_params.f90
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o 00tweedie_params.o

# 1. Calcs_Real Module (Fixes the _evaluaterek_ error)
# Depends on 00tweedie_params.mod because Calcs_Real.f90 uses the tweedie_params_mod.
Calcs_Real.o Calcs_Real.mod: Calcs_Real.f90 00tweedie_params.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o Calcs_Real.o

# Depends on 00tweedie_params.mod because Calcs_Real.f90 uses the tweedie_params_mod.
Calcs_Imag.o Calcs_Imag.mod: Calcs_Imag.f90 00tweedie_params.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o Calcs_Imag.o

# Depends on 00tweedie_params.mod because Calcs_Real.f90 uses the tweedie_params_mod.
Calcs_Solvers.o Calcs_Solvers.mod: Calcs_Solvers.f90 00tweedie_params.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o Calcs_Solvers.o

# Depends on 00tweedie_params.mod because Calcs_Real.f90 uses the tweedie_params_mod.
Calcs_K.o Calcs_K.mod: Calcs_K.f90 Calcs_Solvers.mod 00tweedie_params.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o Calcs_K.o

# 2. Integrands Module (Uses Calcs_Real)
# This is both a compilation rule and a dependency rule.
Integrands.o Integrands_mod.mod: Integrands.f90 Calcs_Real.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o Integrands.o

# 3. gaussianData Module
gaussianData.o gaussian_data_mod.mod: gaussianData.f90 00tweedie_params.mod rprintf_mod.mod
	$(FC) $(FFLAGS) $(FCFLAGS) -fPIC -c $< -o gaussianData.o

# --- Object Dependency Rules ---
# These rules ensure the correct order for object files that consume modules.

# 4. GaussQuadrature.o (Uses Integrands_mod and gaussian_data_mod)
GaussQuadrature.o: Integrands_mod.mod gaussian_data_mod.mod

# 5. TweedieIntegration.o (The main integration routine, depends on many compiled objects)
# It uses Calcs_Real.mod via its dependencies, so we only list the .o files here.
TweedieIntegration.o: Integrands.o accelerate.o GaussQuadrature.o gaussianData.o

# --- Final Linker Flags ---
# If you need to include the RPATH, keep this:
LDFLAGS += -Wl,-rpath,$(R_HOME)/bin/exec/R/lib 

PKG_CPPFLAGS = $(SHLIB_OPENMP_CFLAGS) -I$(R_INCLUDE_DIR)

SDK_PATH := $(shell xcrun --show-sdk-path)
CPPFLAGS += -I$(SDK_PATH)/usr/include/c++/v1   