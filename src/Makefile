LINKFLAGS_FOR=-O2 -march=native -fbounds-check
COMP_FOR = gfortran
OBJS =  vectors.o	\
	main.o		\
        output.o
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c vectors.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -c output.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -c main.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -o R2D2 main.o vectors.o cell.o output.o
clean:;         @rm -f $(OBJS) *.mod *.o
