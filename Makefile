# Executable
EXE = program.exe

# Git Repositary
GIT_URL = git@github.com:mazinmmustafa/MLM_V1.git

# Compilers
CXX = g++
CC = gcc
FC = gfortran

# Compiler Flags
CXXFLG = -Wall -Wextra -pedantic -std=c++17
CFLG = -Wall -Wextra -pedantic -std=c17
F90FLG = -Wall -Wextra -pedantic
F77FLG = -std=legacy

CLIB = -lgfortran -lquadmath -lm

# Optimization
OPT_LEVEL = -O2
OPT = $(OPT_LEVEL) -march=native -mtune=native 
CXXOPT = $(OPT) 
COPT = $(OPT) 
F90OPT = $(OPT) 
F77OPT = -O2 
FLTO = -flto
# FLTO = 

# Directories
BDIR = bin
SDIR = src
ODIR = .obj
DDIR = .dep
HDIR = include
PDIR = data

# Variables
CXXSRC = $(wildcard $(SDIR)/*.cpp)
CXXOBJ = $(patsubst $(SDIR)/%.cpp, $(ODIR)/%.o, $(CXXSRC))
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))
F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))

CXXDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(CXXOBJ))
CDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(COBJ))

CXXDEPFLG = -MMD -MF $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $@) 

# Targets
.DEFAULT_GOAL = build
.PHONY: build run clean clean_all

build: $(BDIR) $(ODIR) $(DDIR) $(PDIR) $(BDIR)/$(EXE) 

run: build
	@echo executing $^:
	@./$(BDIR)/$(EXE)

$(BDIR): 
	@mkdir -pv $@

$(ODIR): 
	@mkdir -pv $@

$(DDIR): 
	@mkdir -pv $@

$(PDIR): 
	@mkdir -pv $@

$(BDIR)/$(EXE): $(CXXOBJ) $(COBJ) $(F90OBJ) $(F77OBJ)
	$(CXX) $(FLTO) -o $@ $^ -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) -g $(CXXFLG) $(CXXOPT) $(CXXDEPFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) -g $(CFLG) $(COPT) $(CDEPFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -g $(F90FLG) $(F90OPT) -o $@ -c $<

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) $(F77OPT) -o $@ -c $<

clean:
	@$(RM) -rv $(BDIR)/$(EXE) $(ODIR)/*.o $(DDIR)/*.d

clean_all: clean 
	@$(RM) -rv $(BDIR) $(ODIR) $(DDIR) .vscode data/*

clean_data:
	@find ./data/ -type f -name '*.pos' | xargs $(RM) -rv
	@find ./data/ -type f -name '*.bin' | xargs $(RM) -rv
	@find ./data/ -type f -name '*.txt' | xargs $(RM) -rv
	@find . -type f -name 'octave-workspace' | xargs $(RM) -rv

.PHONY: valgrind 

valgrind: build
	valgrind --leak-check=full --track-origins=yes ./$(BDIR)/$(EXE) 

.PHONY: git_push git_pull gmsh

git_push: clean_all clean_data
	git add .
	git commit -m 'update' 
	git push --force --set-upstream $(GIT_URL)

git_pull: clean_all
	git pull $(GIT_URL)

gmsh:
	@python3 mesh/generate_mesh.py

# figures:
# 	$(MAKE) -C figures/RCS1/

# Includes
-include $(CXXDEP) $(CDEP)

