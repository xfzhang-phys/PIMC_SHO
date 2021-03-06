# A simple makefile for PIMC code

CC = icpc
CFLAGS = -std=c++11
OPT =  -O3 -ipo -fno-alias -xHost -unroll -no-prec-div
LIBS =
INC = 


EXE = pimc.x
DEP = $(shell find -name "*.h")
SRC = $(shell find -name "*.cpp")
OBJ = $(SRC:%.cpp=%.o)

%.o: %.cpp $(DEP)
	$(CC) $(DFLAGS) $(CFLAGS) $(OPT) -c $< -o $@ $(INC) $(LIBS)

$(EXE): $(OBJ)
	$(CC) $(DFLAGS) $(CFLAGS) $(OPT) -o $(EXE) $(OBJ) $(INC) $(LIBS) 
	rm -rf $(OBJ)

.PHONY: clean clean-all install

clean:
	rm -rf $(OBJ)

clean-all:
	rm -rf $(OBJ) $(EXE)

install:
	cp $(EXE) ../

