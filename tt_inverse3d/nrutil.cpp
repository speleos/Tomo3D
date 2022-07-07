#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

void nrerror(char const* error_text)
/* Numerical Recipes tandard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error (using it was the first one)...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
