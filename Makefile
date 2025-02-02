CC=gcc

q2: 
	$(CC) -o $@ $@.c -lm -lblas -llapacke -llapack   

