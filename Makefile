CXX = g++
CXXFLAGS = -std=c++1z -Ofast -DNODEBUG -W -Wall -Wno-deprecated
LINKFLAGS = -lm

all: NaiveRecompTtoG

NaiveRecompTtoG: NaiveRecompTtoG.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) NaiveRecompTtoG.cpp -o NaiveRecompTtoG

clean:
	rm -f NaiveRecompTtoG *.o *~


