OUT = bobby
BAMTOOLS = /pseudospace1/cumbiej/local
TCLAP = ./tclap

all:
	g++ main.cpp SplatPool.cpp SplatPool.h splats.cpp utils.cpp splats.h utils.h types.h -o $(OUT) -I. -I$(TCLAP) -I$(BAMTOOLS)/include/bamtools -L$(BAMTOOLS)/lib/bamtools -lbamtools -Wall
