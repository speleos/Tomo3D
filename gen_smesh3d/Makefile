# Makefile 
# zeltform.cc is NOT modified to 3D

#CCC = /usr/bin/g++-4.1
#CC  = /usr/bin/gcc-4.1
CCC = g++
CC  = gcc
CFLAGS = -ansi -O3 -pedantic -Wall -pg
HOME = /usr/local/src/tomo3D/gen_smesh3d
PWD = ${HOME}
INCLUDE = -I${PWD}
INSTALL = cp
INSTALLOPT = 

SRC_GEN_SMESH3D_CC = gen_smesh3d.cc zeltform3d.cc \
                     geom3d.cc error.cc util3d.cc
                      
OBJ_GEN_SMESH3D = ${SRC_GEN_SMESH3D_CC:.cc=.o} ${SRC_GEN_SMESH3D_C:.c=.o}

.c.o: 
	$(CCC) $(CFLAGS) -c  $< $(INCLUDE)
.cc.o: 
	$(CCC) $(CFLAGS) -c  $< $(INCLUDE)

gen_smesh3d: $(OBJ_GEN_SMESH3D)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_GEN_SMESH3D) -lm

