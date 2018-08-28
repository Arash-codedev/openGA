# This library is free and distributed under
# Mozilla Public License Version 2.0.


CXX:=g++
#CXX:=x86_64-w64-mingw32-g++
DEBUG_FLAG:= -g -O3
RELESE_FLAG:= -O3 -s -DNDEBUG
CURRENT_FLAGS:= $(RELESE_FLAG)
CURRENT_FLAGS += -std=c++11 -pthread -I./src -Wall -Wconversion -Wfatal-errors -Wextra
BIN=./bin

LIBS:= #-lboost_system -lboost_thread -lboost_chrono
G_LIBS:= -lGL -lGLU -lglut -lGLEW -lSDL -lSDL2main -lSDL2

all:
	@echo "***********************************************"
	@echo Run one of the following commandline examples:
	@echo ""
	@echo make ex_so1
	@echo make ex_so_rastrigin
	@echo make ex_mo1
	@echo make ex_mo_dtlz2
	@echo make ex_iga_colors
	@echo "***********************************************"

ex_so1:
	$(CXX) $(CURRENT_FLAGS) examples/so-1/example_so1.cpp -o $(BIN)/example_so1 $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_so1

ex_so_rastrigin:
	$(CXX) $(CURRENT_FLAGS) examples/so-rastrigin/so-rastrigin.cpp -o $(BIN)/example_so-rastrigin $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_so-rastrigin

ex_mo1:
	$(CXX) $(CURRENT_FLAGS) examples/mo-1/example_mo1.cpp -o $(BIN)/example_mo1 $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_mo1

ex_mo_dtlz2:
	$(CXX) $(CURRENT_FLAGS) examples/mo-dtlz2/mo-dtlz2.cpp -o $(BIN)/example_mo-dtlz2 $(LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/example_mo-dtlz2

ex_iga_colors:
	$(CXX) $(CURRENT_FLAGS) examples/iga-colors/iga-colors.cpp -o $(BIN)/iga-colors $(LIBS) $(G_LIBS)
	@echo "-----------------------------------------------"
	$(BIN)/iga-colors

clean:
	rm ./bin/example_*
