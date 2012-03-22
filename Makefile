OUT = bobby
BAMTOOLS = /Users/abuchanan/bamtools

all:
	g++ main.cpp SplatPool.cpp types.h SplatPool.h -o $(OUT) -I$(BAMTOOLS)/include -L$(BAMTOOLS)/lib -lbamtools
