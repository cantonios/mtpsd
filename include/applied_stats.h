/*

Copyright (c) 2010 C. Antonio Sanchez

This file is part of MTPSD.

MTPSD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MTPSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MTPSD.  If not, see <http://www.gnu.org/licenses/>.

*/

////////////////////////////////////////////////////////////////////
//  NOTE:  Some of this code is a translation of a set of fortran 
//         routines used to calculate the quantiles (ie Percentage 
//         Points) of Gaussian and Chi-squared random variables.
//         The original algorithms are:
//
//          Algorithm AS 91: PPCHI2, DJ Best and DE Roberts (1975)
//          Algorithm AS 32: GAMAIN, GP Bhattacharjee (1970)
//          Algorithm AS 70: GAUINV, RE Odeh and JO Evans (1974)
///////////////////////////////////////////////////////////////////

#define AS_TOL 1e-6    //Accuracy of methods

namespace stats {

// Inverts chi-squared CDF with v degrees of freedom
// to find percentage point x satisfying  P( chi^2_v < x ) = p
// Returns -1 on error
double chisquared_inv(double p, double v);

// Inverts chi-squared CDF for an array (length n) of DOFs v, and 
// places the corresponding percentage points in x 
void chisquared_inv(double p, double *v, unsigned int n, double *x);

// Inverts chi-squared CDF for array of DOFs and PPs
// x filled s.t. x[ii][jj] = x[nv*ii+jj] i.e. ii refers to PP, jj refers to DOF
void chisquared_inv(double *p, double *v, unsigned int np, unsigned int nv, double *x);

// Inverts standard Gaussian CDF to find pp
// x satisfying  P( N < x ) = p
double gauss_inv(double p);

// Computes the incomplete gamma integral:
// I(x,n) = 1/gamma(n) * int ( exp(-t)*t^(n-1), t=0..x)
double gamma_int(double x, double n);

// Computes incomplete gamma integral using series expansion
double gamma_int_series(double x, double n);

// Computes incomplete gamma integral using continued fractions
double gamma_int_fraction(double x, double n);

}
