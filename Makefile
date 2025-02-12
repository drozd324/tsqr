CC = gcc
PACKAGES = -lm -lblas -llapacke -llapack


q2: 
	$(CC) -o $@ $@.c $(PACKAGES)

qr: 
	$(CC) -o $@ $@.c $(PACKAGES)

clean:
	rm qr q2

