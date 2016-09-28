FC = gfortran $(FFLAGS)  
CXX = g++ $(CXXFLAGS) 

p1 = run_IMSRG
p2 = get_HO_energies
p3 = truncate

FFLAGS = -O3 -fopenmp
TFLAGS = -g -O0 -fbounds-check -fopenmp
PFLAGS = -O3 -pg

CXXFLAGS = -std=c++11 -O3

LIBS =  -L/user/local/lib/ -llapack -lblas -larpack -lz

obfiles = bin
modfiles = md

F90SRC=$(wildcard *.f90)
OBJ = $(patsubst %.f90, $(obfiles)/%.o, $(wildcard *.f90))  
OBJ_f = $(patsubst %.f, $(obfiles)/%.o, $(wildcard *.f))  

all: $(OBJ_f) $(OBJ)
	${FC} $^ -o ${p1} -J$(modfiles) ${LIBS}

$(OBJ): | $(obfiles)
$(OBJ_f): | $(obfiles)

$(obfiles):
	@mkdir -p $@
	@mkdir -p $(modfiles)

$(obfiles)/%.o: %.f
	${FC} -c -o $@ $< -J$(modfiles)

$(obfiles)/%.o: %.f90
	${FC} -c -o $@ $< -J$(modfiles)

# nice gift from FEI to detect dependencies automatically
dependencies.mk: $(F90SRC)
	@for f in $^; do \
	    printf "%s:" "$(obfiles)/$${f%.f90}.o"; \
	    awk -v p="$(obfiles)/" \
	        '$$1 == "use" && NF == 2 { printf " %s%s.o",p,$$2 }' "$$f"; \
	    echo; \
	done >$@.tmp; \
	mv $@.tmp $@
 
-include dependencies.mk

clean:
	rm -f ${p1} $
	rm -f ${p2} $
	rm -f ${p3} $
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -rf $(obfiles)
	rm -rf $(modfiles)
	rm -f dependencies.mk

