CPP = g++ -std=c++11
CC = gcc
OPT = -O3 -g -fPIC -Wall -mfpmath=sse -msse4a -fexpensive-optimizations -pipe

CWD 	= $(shell pwd)

SOURCE 	= main.cc


INCLUDE = . ../../matrix/ ../../common/ ../../number_generators/ ../../number_generators/random/ ../../number_generators/low_descrepancy/ ../../statistics/ $(TBB_INCLUDE_PATH)/tbb/

DEFINES = 

OBJECTS = $(SOURCE:.cc=.o)

LIBPATH = $(TBB_LIB_PATH)
LIBRARY = tbb

# All
all:: hjm

# Compile source code to object files
.cc.o :
	$(CPP) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@
        
.c.o :
	$(CC) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@


# Link the object files to build the simulator
hjm : $(OBJECTS) 
	$(CPP) -o $@ $(OBJECTS) $(patsubst %,-L%,$(LIBPATH)) $(patsubst %,-l%,$(LIBRARY))


# Clean targets
clean ::
	$(RM) $(OBJECTS)
	$(RM) hjm
