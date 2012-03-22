OUT = bobby
BAMTOOLS = /Users/abuchanan/bamtools

all:
	g++ main.cpp SplatPool.cpp SplatPool.h splats.cpp utils.cpp splats.h utils.h types.h -o $(OUT) -I$(BAMTOOLS)/include -L$(BAMTOOLS)/lib -lbamtools
