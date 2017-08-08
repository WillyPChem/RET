LIB        = -L. -lfftw3 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = Coupled.x
CXX        = g++

${EXEC}: CoupledDensityMatrix.c 
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} CoupledDensityMatrix.c -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

