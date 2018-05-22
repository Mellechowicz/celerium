CXX=g++
CFLAGS= -std=c++17 -Wall -pedantic -O3 -c 
LDFLAGS=  
#DEBUGFLAG = -D_VERBOSE #uncomment for verbose

BIN = ./obj

MODE = test_potential

SRC =

MKDIR_P = mkdir -p

default: all

obj:
	${MKDIR_P} obj

clean:
	rm -rf $(BIN) *.x data.dat

all: obj	
	@OBJS="";									\
	NAMES="";									\
	for _file in ${SRC};								\
	do										\
		printf "Building: $${_file}.cpp...\t";					\
		$(CXX) $(CFLAGS) $(DEBUGFLAG) $${_file}.cpp -o $(BIN)/$${_file}.o;	\
		echo "done.";								\
		OBJS="$(BIN)/$${_file}.o $${OBJS}";					\
		NAMES="$${_file}.o $${NAMES}";						\
	done;										\
	for _file in ${MODE};								\
	do										\
		printf "Building: $${_file}.cpp...\t";					\
		$(CXX) $(CFLAGS) $(DEBUGFLAG) $${_file}.cpp -o $(BIN)/$${_file}.o;	\
		echo "done.";								\
		printf "Linking: $${NAMES}$${_file}.o\t-> $${_file}.x...\t";		\
		$(CXX) $(LDFLAGS) $${OBJS} $(BIN)/$${_file}.o -o $${_file}.x;		\
		echo "done.";								\
	done

