#include <stdio.h>
#include <math.h>
#include "mtpsd.h"

#ifndef PI
#define PI 3.14159265
#endif

int main( int argc, const char* argv[] ) {

	double f = 2.0;
	double w = 2*PI*f;
	const uint_t n = 200;
	double *x = new double[n];

	for (uint_t i=0; i<n; i++) {
		x[i] = sin(i*w/n);
		printf("s=%f, t=%fs \n", x[i], ((double)i)/n);
	}
	printf("\n");

	mtpsd<double> spectrum(x, n, 1.5);
	try {
		spectrum.compute();
	} catch (ERR e) {
		printf("ERROR: %s\n", e.getmsg());
	}

	for (uint_t i=0; i<n; i++) {
		printf("S=%f, f=%fHz \n", spectrum(i), spectrum.freq(i));
	}

}
