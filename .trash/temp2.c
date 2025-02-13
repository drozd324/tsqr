# include <stdio.h>

struct example {
    double* c;
	int m;
};

int main() {
    struct example e;
    printf("&e = %p\n", &e);
    printf("&e.c = %p\n", &(e.c));
    printf("&(e.c) = %p\n", &(e.c));
    printf("&(e.m) = %p\n", &(e.m));
    //printf("Size of example structure: %ld bytes\n", sizeof(e));
    return 0;
}
