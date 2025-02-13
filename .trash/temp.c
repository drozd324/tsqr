#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main( void )
{
    int a;

    a = fmin( 1, 0.10 );
    printf( "The value is: %d\n", a );
    return EXIT_SUCCESS;
}
