#
# compiler
#
CC        = g++
#CC       	  = clang++
CXX_FLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG


TARGET = test

all: $(TARGET)

%.o:%.cc
	$(CC) -c $(CXX_FLAGS) $(CUDD_INCa) $(CUDD_INCb) $< -o $@

$(TARGET): $(TARGET).o
	$(CC) $(CXX_FLAGS) -o $(TARGET) $(TARGET).o $(CUDD_LIB_PATH) $(CUDD_LIBS)


clean:
	rm -f ./$(TARGET)  ./$(TARGET).o
