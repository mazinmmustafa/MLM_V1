# Executable
EXE = program.exe

# Git Repositary
GIT_URL = git@github.com:mazinmmustafa/MLM_V1.git master

# Compiler
CC = gcc
FC = gfortran

# Compiler Library
CLIB = -lgfortran -lquadmath -lm

# Compiler Flags
CFLG = -Wall -Wextra -std=c11 -pedantic
OPT_LVL = -O2
OPT = $(OPT_LVL) -march=native -mtune=native

# Directories
BDIR = bin
HDIR = include
SDIR = src
ODIR = .obj
DDIR = .dep

# Variables
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
CDEP = $(patsubst $(SDIR)/%.c, $(DDIR)/%.d, $(CSRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))
F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))

# Dependency Flag
DFLG = -MMD -MF $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $@)

# Targets
.PHONY: run build

build: $(ODIR) $(DDIR) $(BDIR) $(BDIR)/$(EXE)

$(ODIR):
	@mkdir -pv $@

$(DDIR):
	@mkdir -pv $@

$(BDIR):
	@mkdir -pv $@

run: build 
	@./$(BDIR)/$(EXE)

$(BDIR)/$(EXE): $(F77OBJ) $(F90OBJ) $(COBJ)
	$(CC) -g -o $@ $^ -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) -g $(DFLG) $(OPT) $(CFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(OPT) -o $@ -c $< 

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(OPT) -o $@ -c $< 

-include $(CDEP)

# Debug
.PHONY: debug
debug:
	@echo $(SDIR)
	@echo $(CSRC)
	@echo $(COBJ)

# Clean
.PHONY: clean
clean:
	@$(RM) -rv $(EXE) $(COBJ) $(CDEP)

clean_all: clean
	@$(RM) -rv $(BDIR) $(ODIR) $(DDIR)

.PHONY: valgrind 

valgrind: build
	valgrind --leak-check=full --track-origins=yes ./$(BDIR)/$(EXE) 

.PHONY: git_push git_pull

git_push: clean clean_all
	git add .
	git commit -m 'update' 
	git push --force --set-upstream $(GIT_URL)

git_pull: clean clean_all
	git pull $(GIT_URL)