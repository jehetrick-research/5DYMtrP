#
#
# -DLINUX if on Linux/Unix NO -DLINUX on cygwin
# -DUSE_DOUBLE_PRECISION for double precision (the most efficient option)
# -DSSE2 -O3 for pentium4 optimizations
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
# NOTE: I hate makefiles and I do not use them. 
# If you know how to write a better makefile send it to me please!
#

#CC     = g++ -I../Libraries 
CC     = mpic++ -I../Libraries 
# Basic compilation flags (-O2 for speed)
#CFLAGS = -DLINUX -O3
# Flags for double precision and sse (P4 optiomizations)
#CFLAGS = -DLINUX -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CFLAGS = -DPARALLEL -DLINUX -DUSE_DOUBLE_PRECISION -DSSE2 -O2
# -Wno-write-strings suppresses depricated "conversion from string constant to ‘char*’" warning
CFLAGS = -DPARALLEL -DLINUX -O2 -Wno-write-strings


pureSUN::
	${CC} ${CFLAGS} pureSUN.cpp -o pureSUN


su3metrop:: su3metrop.cpp metropolis.cpp
	${CC} ${CFLAGS} su3metrop.cpp -o su3metrop

su3metropP:: su3metropP.cpp ploop.cpp metropolis.cpp
	${CC} ${CFLAGS} su3metropP.cpp -o su3metropP

su3metropPls:: su3metropPls.cpp ploop.cpp metropolis.cpp
	${CC} ${CFLAGS} su3metropPls.cpp -o su3metropPls

tests::
	${CC} ${CFLAGS} tests.cpp -o tests

example02::
	${CC} ${CFLAGS} example02.cpp -o example02

example03::
	${CC} ${CFLAGS} example03.cpp -o example03

su3trP:: su3trP.cpp ploop3.cpp
	${CC} ${CFLAGS} su3trP.cpp -o su3trP

su3dev::
	${CC} ${CFLAGS} su3dev.cpp -o su3dev

su3dev2:: su3dev2.cpp ploop2.cpp
	${CC} ${CFLAGS} su3dev2.cpp -o su3dev2

su3dev3:: su3dev3.cpp ploop2.cpp
	${CC} ${CFLAGS} su3dev3.cpp -o su3dev3

su3dev4:: su3dev4.cpp ploop3.cpp
	${CC} ${CFLAGS} su3dev4.cpp -o su3dev4

su3devtest:: su3devtest.cpp ploop3.cpp
	${CC} ${CFLAGS} su3devtest.cpp -o su3devtest

