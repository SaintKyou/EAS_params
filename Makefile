CC=g++
CFLAGS=-c 

LDFLAGS= -lpthread -pthread -fcommon -std=c++11
SOURCES=main.cpp metropolis.cpp matrix.cpp functions.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.o main
