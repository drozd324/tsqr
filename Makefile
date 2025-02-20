CC    = gcc -O0
MPICC = mpicc -O0
DEBUG = -g -fsanitize=address -Wall -Wextra -lefence #$(" ") #
LIBS  = -lm -lblas -llapacke -llapack  

all: q2 q3

q2: q2.c mat_tools.o tsqr.o
	$(MPICC) -o q2 q2.c mat_tools.o tsqr.o $(LIBS) $(DEBUG) 

q3: q3.c mat_tools.o tsqr.o
	$(MPICC) -o q3 q3.c mat_tools.o tsqr.o $(LIBS) $(DEBUG)

mat_tools.o: mat_tools.c mat_tools.h
	$(CC) -c mat_tools.c $(LIBS) $(DEBUG)

tsqr.o: tsqr.c tsqr.h mat_tools.c mat_tools.h
	$(MPICC) -c tsqr.c $(LIBS) $(DEBUG)

clean:
	rm *.o q2 q3
