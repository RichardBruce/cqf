CPP = g++
CC = gcc
CC_OPT = -O3 -g -std=c++0x -fPIC -Wall -mfpmath=sse -msse4a -fexpensive-optimizations -pipe
LD_OPT = -static

CWD 	= $(shell pwd)

SOURCE 	= main.cc

INCLUDE = . ../number_generators/random ../common

DEFINES = 

OBJECTS = $(SOURCE:.cc=.o)

LIBPATH = 
LIBRARY = 

# All
all:: analytic

# Compile source code to object files
.cc.o :
	$(CPP) $(CC_OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@
        
.c.o :
	$(CC) $(CC_OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@


# Link the object files to build the simulator
analytic : $(OBJECTS) 
	$(CPP) -o $@ $(LD_OPT) $(OBJECTS) $(patsubst %,-L%,$(LIBPATH)) $(patsubst %,-l%,$(LIBRARY))


# Clean targets
clean ::
	$(RM) $(OBJECTS)
	$(RM) analytic
