OUT = bobby
BAMTOOLS_INC = /pseudospace1/cumbiej/local/include/bamtools
BAMTOOLS_LIB = /pseudospace1/cumbiej/local/lib/bamtools

all:
	g++ bobby.cpp -o $(OUT) -I $(BAMTOOLS_INC) -L$(BAMTOOLS_LIB) -lbamtools
