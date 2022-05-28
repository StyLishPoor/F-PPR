CC      = g++
CPPFLAGS= -Wno-deprecated -O3 -c -m64 -fopenmp -march=native -std=c++17 -DGCC -DRelease
#CPPFLAGS= -Wall -O3 -c -m64 -fopenmp -march=native -std=c++17 -DGCC -DRelease
#LDFLAGS = -static -O3 -m64 -fopenmp
LDFLAGS = -O3 -m64 -fopenmp
SOURCES = main-dev-pull.cpp Util.cpp Graph.cpp UnitHeap.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE=Frontier-dev-pull

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o

