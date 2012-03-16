OUT = bobby
BAMTOOLS = 

all:
	g++ bobby.cpp -o $(OUT) -I$(BAMTOOLS)/include -L$(BAMTOOLS)/lib -lbamtools
