CPP = g++
CC = gcc
CC_OPT = -O3 -g -std=c++0x -fPIC -Wall -mfpmath=sse -msse4a -fexpensive-optimizations -pipe
LD_OPT = 

CWD 	= $(shell pwd)

SOURCE 	= main.cc

INCLUDE = . ../matrix/ ../common/ ${BOOST_INCLUDE_PATH} ${TBB_INCLUDE_PATH}

DEFINES = 

OBJECTS = $(SOURCE:.cc=.o)

LIBPATH = ${BOOST_LIB_PATH} ${TBB_LIB_PATH}
LIBRARY = tbb boost_thread

# All
all:: pde_solver

# Compile source code to object files
.cc.o :
	$(CPP) $(CC_OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@
        
.c.o :
	$(CC) $(CC_OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@


# Link the object files to build the simulator
pde_solver : $(OBJECTS) 
	$(CPP) -o $@ $(LD_OPT) $(OBJECTS) $(patsubst %,-L%,$(LIBPATH)) $(patsubst %,-l%,$(LIBRARY))


# Clean targets
clean ::
	$(RM) $(OBJECTS)
	$(RM) pde_solver
