CXX = g++
CFLAGS = -std=c++17 -O3 -Wall -pedantic -c 
LDFLAGS = -lm

BIN = ./obj

MAIN_PREFIX = main
MODE = vector

MKDIR_P = mkdir -p

$(MODE): obj $(BIN)/$(MAIN_PREFIX).o
	$(CXX) $(BIN)/$(MAIN_PREFIX).o $(LDFLAGS) -o $(MODE)

$(BIN)/$(MAIN_PREFIX).o: $(MAIN_PREFIX).cpp
	$(CXX) $(CFLAGS) $(MAIN_PREFIX).cpp -o $(BIN)/$(MAIN_PREFIX).o 

obj:
	${MKDIR_P} obj

clean:
	rm -rf $(MODE) $(BIN)

run:
	./$(MODE)

all:	clean obj $(MODE) run
