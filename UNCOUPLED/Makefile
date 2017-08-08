LIB        = -L. -lfftw3 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = NP.x
CXX        = g++

${EXEC}: NP.c 
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} NP.c -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

