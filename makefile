FC      := gfortran
EXE     := octopus

OBJDIR  := build/obj
MODDIR  := build/mod

# ---- flags ----
FFLAGS  := -O3 -Wall -Wextra -fimplicit-none -fopenmp -ffree-line-length-none
MODFLAGS := -J$(MODDIR) -I$(MODDIR)
LDFLAGS := -fopenmp

SRC     := parameters.f90 dynamics.f90 emission.f90 functions.f90 \
           initial_condition.f90 mainprogram.f90 method.f90 \
           model.f90 operations.f90

OBJS    := $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))

.PHONY: all clean dirs
all: $(EXE)

dirs:
	mkdir -p $(OBJDIR) $(MODDIR)

$(EXE): dirs $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

$(OBJDIR)/%.o: %.f90 | dirs
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

clean:
	rm -rf build $(EXE)