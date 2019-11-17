# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -O3 lib.cpp
MPICPPflags = mpic++ -std=c++11 -O3 lib.cpp
# Here we define the library functions we nee
# LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable

MPIPROG = mpi_ising
PROG = ising
TESTPROG = testcode
all: ${MPIPROG} ${PROG} ${TESTPROG}

${PROG}      : main.o isingmodel.o
			    ${CPPflags} main.o isingmodel.o -o ${PROG}
main.o       : main.cpp
				${CPPflags} -c main.cpp
isingmodel.o : isingmodel.cpp
				${MPICPPflags} -c isingmodel.cpp

${MPIPROG} : mainMPI.o isingmodel.o
				${MPICPPflags} mainMPI.o isingmodel.o -o ${MPIPROG}
mainMPI.o    : mainMPI.cpp
				${MPICPPflags} -c mainMPI.cpp

${TESTPROG} :	 	tests-main.o  test_functions.o isingmodel.o
			${CPPflags} tests-main.o test_functions.o isingmodel.o ${LIB} -o ${TESTPROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test_functions.cpp
		        	${CPPflags} -c test_functions.cpp

jacobi.o :		isingmodel.cpp
		        	${CPPflags} -c isingmodel.cpp

