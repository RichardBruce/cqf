CPP = g++
OPT = -O3 -g -fPIC -Wall -mfpmath=sse -msse4a -fexpensive-optimizations -pipe

INCLUDE = ../number_generators/random/ ../common/
DEFINES = 
SOURCE = monte_carlo.cc
LIBPATH = 
LIBRARY = 

OBJECTS = $(SOURCE:.cc=.o)

# Compile source code to object files
.cc.o :
	$(CPP) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@
        
.c.o :
	$(CC) $(OPT) $(patsubst %,-D%,$(DEFINES)) $(patsubst %,-I%,$(INCLUDE)) -c $< -o $@

# MC make file
mc : $(OBJECTS) 
	$(CPP) -o $@ $(OBJECTS) $(patsubst %,-L%,$(LIBPATH)) $(patsubst %,-l%,$(LIBRARY))

clean ::
	$(RM) $(OBJECTS)
	$(RM) mc
