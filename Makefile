SERVER=FX

###############################################
ifeq ($(SERVER), LOCAL)
FC = mpif90
FFLAGS = -Wall -Wno-tabs
endif
###############################################
ifeq ($(SERVER), FX)
FC = mpifrtpx
FFLAGS = -Kfast,openmp,parallel -fPIC -Nalloc_assign
endif


UPDATE = update
SRC_FILES := $(wildcard ./*.F90)
OBJ_FILES := $(notdir $(SRC_FILES:.F90=.o))
SRC_DIR   := $(dir $(SRC_FILES))
VPATH := $(SRC_DIR)

$(UPDATE) : $(OBJ_FILES)
	$(FC) $(FFLAGS) -o a.out $(OBJ_FILES) && \
	touch $(UPDATE)

%.o: %.F90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -rf $(OBJ_FILES) a.out $(UPDATE) *.mod

main.o: param_def.o change_judge.o
upgrade_param_set.o: param_def.o
