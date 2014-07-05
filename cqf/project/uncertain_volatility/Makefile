CPP = g++ -std=c++11
CC = gcc
OPT = -O3 -g -fPIC -Wall -mfpmath=sse -msse4a -fexpensive-optimizations -pipe

CWD 	= $(shell pwd)

SOURCE 	= main.cc fd_solver.cc downhill_simplex.cc quasi_newton_solver.cc metropolis.cc


INCLUDE = . ../../matrix/ ../../common/ ../../number_generators/ ../../number_generators/random/

DEFINES = 

OBJECTS = $(SOURCE:.cc=.o)

LIBPATH = 
LIBRARY = 

# All
all:: uncertain_volatility

# Compile source code to object files
.cc.o :
	$(CPP) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@
        
.c.o :
	$(CC) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@


# Link the object files to build the simulator
uncertain_volatility : $(OBJECTS) 
	$(CPP) -o $@ $(OBJECTS) $(patsubst %,-L%,$(LIBPATH)) $(patsubst %,-l%,$(LIBRARY))


# Clean targets
clean ::
	$(RM) $(OBJECTS)
	$(RM) uncertain_volatility
        
