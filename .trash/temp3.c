# include <stdio.h>

struct example {
    double* c;
	int m;
};

void func(struct example e){
	printf("%p\n", &e);
	//printf("%d\n", e->m);
	printf("%d\n", e.m);
}

int main() {
    struct example e;
	//(*e).m = 69;
	//e->m = 69;
	e.m = 69;
		
	//int a = 1;
	//int* b = 8;
	//int c[3] = {1, 2, 3};
	
	func(e);

    return 0;
}
