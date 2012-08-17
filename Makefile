#ViennaRNA package location
ViennaRNA=~/software/include/ViennaRNA
ViennaLIB=~/software/lib/libRNA.a

CPP = g++
CC = gcc
CPPFLAGS =  -O2 -g -Wall -std=c++0x -fexceptions -Wno-write-strings
CFLAGS =  -O2 -g -Wall -fexceptions -Wno-write-strings
LFLAGS = -O2 -fopenmp

OBJ = RNAlandmap_cmdline.o\
			main.o\
			move_set.o\
			RNAlandmap.o

DIRS = -I $(ViennaRNA)

LIBS = $(ViennaLIB)

all: $(OBJ)
	$(CPP) $(LFLAGS) $(DIRS) $(OBJ) $(LIBS) -o RNAlocmin
	rm -f $(OBJ)

RNAlandmap_cmdline.h RNAlandmap_cmdline.c: RNAlandmap.ggo
	gengetopt -i RNAlocmin.ggo

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(DIRS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(DIRS) -c $<

clean:
	rm -f $(OBJ)
	rm -f RNAlocmin
	rm -f RNAlocmin_cmdline.c RNAlocmin_cmdline.h

