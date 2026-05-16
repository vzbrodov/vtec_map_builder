FC ?= gfortran
BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
FFLAGS ?= -O2 -Wall -Wextra -std=f2008 -Isrc -I$(MOD_DIR) -J$(MOD_DIR)

SRC = src/my_prec.f90 \
      src/string_utils.f90 \
      src/station_data_mod.f90 \
      src/spherical_harmonics_mod.f90 \
      src/least_squares_mod.f90 \
      src/data_io_mod.f90 \
      src/vtec_model_mod.f90 \
      src/main.f90

OBJ = $(patsubst src/%.f90,$(OBJ_DIR)/%.o,$(SRC))
BIN = build_vtec_maps

.PHONY: all clean run

all: $(BIN)

$(BIN): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

$(OBJ_DIR)/%.o: src/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

run: $(BIN)
	./$(BIN)

clean:
	rm -rf $(BIN) $(BUILD_DIR) src/*.o src/*.mod *.mod
