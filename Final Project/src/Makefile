 
native: CC=g++
native: CFLAGS = -std=c++14 -Wall -pedantic -Werror -D_SEQ_
native: OUT=LBSim_SEQ

parallel: CC=mpic++
parallel: CFLAGS = -std=c++14 -Wall -pedantic -Werror -D_MPI_
parallel: OUT=LBSim_MPI

.cpp.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
		visual.o\
		readconfig.o\
		simulation.o\
		main.o

native: $(OBJ)
	$(CC) $(CFLAGS) -o $(OUT) $(OBJ)  -lm

parallel: $(OBJ)
	$(CC) $(CFLAGS) -o $(OUT) $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ)


helper.o      	: helper.h 
visual.o      	: helper.h
readconfig.o 	: readconfig.h LB_consts.h helper.h 
simulation.o 	: simulation.h LB_consts.h visual.h
main.o        	: simulation.h readconfig.h