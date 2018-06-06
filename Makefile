CXX=g++
CXXFLAGS=-g -Wall -std=c++11
PROG=polytest
SRC=PolyFit.cpp test.cpp

all: $(PROG)

$(PROG) : $(SRC)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(SRC)

clean : $(PROG)
	rm $(PROG)