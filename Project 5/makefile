# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -O3 lib.cpp
#MPICPPflags = mpic++ -O3 lib.cpp
# Here we define the library functions we nee
 LIB = -larmadillo 
# Here we define the name of the executable

#MPIPROG = mpi_
PROG = PDEsolver
TESTPROG = testcode
all: ${PROG} ${TESTPROG}

${PROG}      : main.o PDEsolver.o
			    ${CPPflags} main.o PDEsolver.o ${LIB}-o ${PROG}


${TESTPROG} :	 	tests-main.o  test_functions.o PDEsolver.o
			${CPPflags} tests-main.o test_functions.o PDEsolver.o ${LIB} -o ${TESTPROG}

main.o       : main.cpp
				${CPPflags} -c main.cpp
PDEsolver.o : PDEsolver.cpp
				${CPPflags} -c PDEsolver.cpp

tests-main.o :	tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test_functions.cpp
		        	${CPPflags} -c test_functions.cpp

# ${MPIPROG} : mainMPI.o isingmodel.o
# 				${MPICPPflags} mainMPI.o isingmodel.o -o ${MPIPROG}
# mainMPI.o    : mainMPI.cpp
# 				${MPICPPflags} -c mainMPI.cpp
