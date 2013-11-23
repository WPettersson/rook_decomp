CPP=g++
CC_FLAGS=-O2 -Wall
#CC_FLAGS=-g -Wall
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXEC=decomp
%.o: %.cpp
	$(CPP) -c $(CC_FLAGS) $< -o $@

$(EXEC): $(OBJECTS)
	$(CPP) $(CC_FLAGS) $(OBJECTS) -o $(EXEC)

clean: 
	rm $(OBJECTS) $(EXEC)
