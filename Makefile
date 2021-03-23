FC = mpif90
UPDATE = update
FFLAGS = -Wall -Wno-tabs
SRC_FILES := $(wildcard ./*.F90)
OBJ_FILES := $(notdir $(SRC_FILES:.F90=.o))
SRC_DIR   := $(dir $(SRC_FILES))
VPATH := $(SRC_DIR)

$(UPDATE) : $(OBJ_FILES)
	$(FC) -o a.out $(OBJ_FILES) && \
	touch $(UPDATE)

%.o: %.F90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -rf $(OBJ_FILES) a.out $(UPDATE) *.mod

main.o: param_def.o change_judge.o
upgrade_param_set.o: param_def.o
