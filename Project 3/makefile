# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3 lib.cpp gauss-laguerre.cpp
# Here we define the library functions we nee
# LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= testcode
${PROG} :	   	tests-main.o  test_functions.o Gauss.o
			${CPPflags} tests-main.o test_functions.o Gauss.o ${LIB} -o ${PROG}

tests-main.o :			tests-main.cpp
		        	${CPPflags} -c tests-main.cpp

test-functitons.o :			test_functions.cpp
		        	${CPPflags} -c test_functions.cpp

jacobi.o :		Gauss.cpp
		        	${CPPflags} -c Gauss.cpp