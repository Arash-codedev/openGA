# This library is free and distributed under
# Mozilla Public License Version 2.0.


CXX:=g++
DEBUG_FLAG:= -g -O3
RELESE_FLAG:= -O3 -s -DNDEBUG -DARMA_NO_DEBUG
CURRENT_FLAGS:= $(RELESE_FLAG)
CURRENT_FLAGS += -std=c++11 -pthread -I./core
BIN=./bin

LIBS:= -lboost_system -lboost_thread -lboost_chrono
LIBS+= -larmadillo

all:
	@echo "***********************************************"
	@echo Run one of the following commandline examples:
	@echo ""
	@echo make ex_so1
	@echo make ex_mo1
	@echo make ex_iga1
	@echo "***********************************************"

ex_so1:
	$(CXX) $(CURRENT_FLAGS) examples/example_so1/example_so1.cpp -o $(BIN)/example_so1 $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_so1

ex_mo1:
	$(CXX) $(CURRENT_FLAGS) examples/example_mo1/example_mo1.cpp -o $(BIN)/example_mo1 $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_mo1

ex_iga1:
	$(CXX) $(CURRENT_FLAGS) examples/example_iga1/example_iga1.cpp -o $(BIN)/example_iga1 $(LIBS) -lGL -lGLU -lglut -lGLEW -lSDL -lSDL2main -lSDL2
	@echo "-----------------------------------------------"
	$(BIN)/example_iga1

clean:
	rm ./bin/example_*
