CPP       = /opt/nvidia/hpc_sdk/Linux_x86_64/21.3/compilers/bin/nvc++
OBJ     = o
EXE     = out
RUN     =
CCFLAGS  = -fast
ACCFLAGS = -Minfo -acc $(OPT)
INC =
#INC += -I/home/nikita.mushak/apps/apps
UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
OBJ     = obj
EXE     = exe
endif

all: build run verify

build: binModel.cpp binomialOPENACC.cpp
	$(CPP) $(INC) -c $(CCFLAGS) binomialOPENACC.cpp
	$(CPP) $(INC) $(CCFLAGS) $(ACCFLAGS) -o binModel.$(EXE) binModel.cpp binomialOPENACC.$(OBJ)

run: binModel.$(EXE)
	$(RUN) ./binModel.$(EXE)

verify:


clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.$(OBJ) *.dwf *.pdb prof


###############################################################################
