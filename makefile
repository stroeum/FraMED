#!bin/bash

.PHONY: clean
.PHONY: print

# Setup directories & folders (ADIR: absolute directory, CDIR: current directory)
ADIR = $(PWD)
BDIR = bin
CDIR = $(notdir $(ADIR))
IDIR = include
LDIR = lib
ODIR = build
RDIR = results
SDIR = src

# Setup files & libs
DEPS = $(wildcard $(IDIR)/*.h)
SRCS = $(wildcard $(SDIR)/*.cpp)
OBJS = $(subst $(SDIR), $(ODIR), $(SRCS:.cpp=.o)) 
LIBS = -lm
EXE  = FraMED

# Setup compiler
CXX = g++
CXXFLAGS = -Wall -Os -std=c++0x -I$(IDIR)

# Setup the q
MAIL    = #rioussej@erau.edu
RUNTIME = 04:23:59:59
NODES   = 1
PPN     = 1
PROCS   = $(shell echo $(NODES)*$(PPN) | bc)
MEMORY  = 8gb
JOBID   = $(shell cat $(BDIR)/JOB.ID | tail -n+2)

# Compile
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	mkdir -p $(ODIR)	
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

# Run
run:
	mkdir -p $(RDIR)
	./$(EXE)

qrun:
	mkdir -p $(RDIR)
	msub \
	  -l walltime=$(RUNTIME),nodes=$(NODES):ppn=$(PPN),mem=$(MEMORY)\
	  -m abe -M $(MAIL)\
	  -o $(ADIR)/$(BDIR)/log_$(CDIR).out\
	  -e $(ADIR)/$(BDIR)/log_$(CDIR).err\
	  -v EXE=$(EXE),ADIR=$(ADIR),BDIR=$(BDIR),CDIR=$(CDIR),PROCS=$(PROCS)\
	  start.job > $(ADIR)/$(BDIR)/JOB.ID 

prun: 
	mpirun -np $(PROCS) ./$(EXE)

# Stop
qstop: 
	qdel $(JOBID)

# Clean
clean:
	rm -rfv $(ODIR)/*.o *~ $(IDIR)/*~

reset: 
	rm -rfv \
	  $(EXE) \
	  $(BDIR)/JOB.ID $(BDIR)/*log* \
	  $(RDIR)/*.?at $(RDIR)/*.avi $(RDIR)/*.eps $(RDIR)/*.ps $(RDIR)/*.pdf $(RDIR)/*.txt

# Display
print:
	$(CDIR)
