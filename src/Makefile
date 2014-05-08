include $(ALPS_ROOT)/share/alps/include.mk
include compiler.flags

all: main

OBJ = io.o interaction_expansion.o splines.o fastupdate.o model.o solver.o worm.o main.o 

main: $(OBJ) observables.o
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(RMFLAGS) *.o

