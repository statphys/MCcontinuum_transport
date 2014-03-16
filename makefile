CC=g++ 
CFLAGS=-c 
LDFLAGS=
SOURCES=mymain.cpp lattice3dMC.cpp conductivity.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=provamc

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean: rm -f *.o core

