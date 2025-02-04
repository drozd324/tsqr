CC=gcc

q2: 
	$(CC) -o $@ $@.c -lm -lblas -llapacke -llapack

qr: 
	$(CC) -o $@ $@.c -lm -lblas -llapacke -llapack

clean:
	rm qr q2

