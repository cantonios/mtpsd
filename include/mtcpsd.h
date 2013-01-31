/*

Copyright (c) 2013 C. Antonio Sanchez

This file is part of mtcpsd.

mtcpsd is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mtcpsd is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mtcpsd.  If not, see <http://www.gnu.org/licenses/>.

*/

/////////////////////////////////////////////////////////////////////////////////////////////
//  mtcpsd.cpp contains all functions required for calculating the multitaper
//            spectrum estimate for a one-dimensional time-series (real or complex)
//            
//            An interface class is provided ( class mtcpsd{} ) to simplify use.
//            The individual routines are also supplied that work directly on 
//            arrays if you wish to have more control over memory use.
//
/////////////////////////////////////////////////////////////////////////////////////////////

#define mtcpsd_VERSION 1.0

#include "dpss.h"           // methods for computing dpss (so user doesn't have to include it manually)
#include "applied_stats.h"  // chi_squared statistics for confidence intervals

#include <complex.h>
#ifndef complex
#define complex             // in msys, complex.h doesn't define "complex", needed by fftw
#endif
#include <fftw3.h>
#include <cmath>        

#define ZERO_TOL 1e-14      // effectively zero
#define ADAPT_TOL 1e-7      // tolerance for adaptive weights

#if !defined(UINT_T)
#define UINT_T
typedef unsigned int uint_t;   //unsigned integer
#endif

enum WEIGHT_METHOD {EIGEN, EQUAL};
enum DATA_TYPE {REAL_DATA, COMPLEX_DATA};
enum CONF_BOUND {LOWER, UPPER};

// Constants involved in computing psd
struct mtcpsd_workspace{
    uint_t n;                         // length of data sequence
    double nW;                        // time-bandwidth product
    WEIGHT_METHOD weight_method;      // specifies type of weighting
    uint_t N;                         // length of Fourier Transform
    uint_t K;                         // number of tapers (eigenspectra)
    double Fs;                        // sampling Frequency
    bool remove_mean;                 // specifies whether to remove the mean in calculations
    
    DATA_TYPE dtype;                  // specifies whether data is real or complex
    uint_t nwk;                       // number of independent weights per eigenspectrum (usually 1 or N)
    
    // Sets default values
    mtcpsd_workspace(): n(0), nW(2), weight_method(EIGEN), N(0), K(3), Fs(1), remove_mean(true), dtype(REAL_DATA), nwk(0){ }
    mtcpsd_workspace(uint_t n_, double nW_):n(n_),nW(nW_), weight_method(EIGEN), N(n_),K(floor(2*nW_)-1), Fs(1), remove_mean(true), dtype(REAL_DATA), nwk(n_){ }
};

void fix_workspace(mtcpsd_workspace &work);      //fixes possibly invalid entries

template <class T>
class mtcpsd{
    
public:    
    mtcpsd(T *x, T *y, uint_t n, double nW);     //uses dpss{} class to compute tapers
    mtcpsd(T *x, T *y, mtcpsd_workspace work);
    mtcpsd(T *x, T *y, const double *tapers, const double *lambda, uint_t n, uint_t K);      //allows k custom tapers of length n
    mtcpsd(T *x, T *y, const double *tapers, const double *lambda, mtcpsd_workspace work);
          
    ~mtcpsd();
    
    void compute();                          // Computes Spectrum with current information
    
    fftw_complex operator()(uint_t ii);            // returns the ith value of the spectrum estimate
    fftw_complex operator()(uint_t k, uint_t ii);  // returns the ith value of eigenspectrum k
    fftw_complex eig_coeff_x(uint_t k, uint_t ii);   //returns the ith value of eigencoefficient k (Jkx)
    fftw_complex eig_coeff_y(uint_t k, uint_t ii);   //returns the ith value of eigencoefficient k (Jky)
    
    double freq(uint_t ii);                 // returns ith frequency
    
    double dof(uint_t ii);                  // returns equivalent degrees of freedom of estimate at frequency ii
    double wt(uint_t kk, uint_t ii);        // returns w_kk(ii), the kkth weight at frequency ii
    double conf_int(uint_t ii, CONF_BOUND side, double p);     //returns either lower or upper bound of px100% confidence interval at frequency ii
    double conf_factor(uint_t ii, CONF_BOUND side, double p);  //returns either lower or upper confidence factor (conf_int=S[ii]*conf_factor)
    
    double lambda(uint_t kk);                //returns the kth energy concentration, seql<= k <= sequ
    double taper(uint_t kk, uint_t ii);      //returns the iith value of sequence seql<= k <= sequ
    
    uint_t length();                        // returns N, length of S
    uint_t size(int dim);                   //if dim < 1 returns length of S, otherwise returns K
    const fftw_complex *pS();               //returns pointer to Spectrum S
    
    mtcpsd_workspace getinfo();              //returns information about this mtcpsd class
      
protected:
    
    void init(T *x, T *y, uint_t n, uint_t K);     //initializes variables
    void init(T *x, T *y, mtcpsd_workspace work);
    void create_arrays();                       //allocates memory
    
    bool tapers_ready;      // true if tapers are already computed or supplied, false otherwise
    
    T *x;                   // Original Data
    T *y;					// Original Data
    fftw_complex *S;        // Final spectrum
    fftw_complex *Jkx;      // Eigencoefficients, Jkx'*Jky = Sk
    fftw_complex *Jky;
    double *wk;             // weights (unity, energy concentrations);
    double *h;              // data tapers
    double *l;              // taper energies
    mtcpsd_workspace info;  // holds information about the PSD
    
};

// EIGENSPECTRA
// computes the individual eigenspectra of data x (double or fftw_complex), and puts in Sk.  Jk is the set of
// eigencoefficients s.t. |Jk|^2 = Sk.
// n is the length of the data
// N is the length of the FFT
// tapers are the data tapers, with lambda energy concentrations in [-W W] (assumes dpss)
// K is the number of tapers
// remove_mean specifies whether to remove sample mean before calculations
// All arrays in row-major format, where rows correspond to the kth taper/eigenspectra
template <class T>
void eigenspecta(T *x, T *y, uint_t n, const double *h, const double *l, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx, fftw_complex *Jky, fftw_complex *Sk);
template <class T>
void eigencoeffs(T *x, uint_t n, const double *h, const double *l, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx);

// COMBINE_EIGSPEC
// combine eigenspectra using weights wk, with size (K x nwk)
// nwk is the # of unique weights per eigenspectra, is usually nwk = 1 (same weight across all freq) or N (different across all f)
// Weights must be normalized so they sum to one at each frequency
void combine_eig_spec( const fftw_complex *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S);
void combine_eig_coeffs( const fftw_complex *Jkx, const fftw_complex *Jky, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S);

// Combines eigenspectra, fills S, but also computes the normalized absolute difference from
// the previous.  This difference = || S_new - S_old ||/N is returned in diff.
// NOTE: S is both an input AND output in this procedure.
void combine_eig_spec( const fftw_complex *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S, double &diff);
void combine_eig_coeffs( const fftw_complex *Jkx, const fftw_complex *Jky, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S, double &diff);

// DEGREES_OF_FREEDOM
// Calculates the equivalent degrees of freedom of the estimate (assuming large n)
// such that S_estimate ~ S_true chi^2_v / v,
double degrees_of_freedom(uint_t ii, double *wk, uint_t nwk, uint_t K);     // at frequency ii
void degrees_of_freedom( double *wk, uint_t nwk, uint_t K, double *v);      // fills vector v

// CONFIDENCE_INTERVAL
// Computes the px100% confidence interval of the spectrum estimate S
// v is the degrees of freedom (vector of length nv)  (nv = 1 for constant dof across specturm, nv=N for varying)
// Sc has length 2*N, and stores the lower then upper confidence bounds in row-major format.
void confidence_interval( double p, fftw_complex *S, uint_t N, double *v, uint_t nv, double *Sc);
void confidence_interval( double p, fftw_complex *S, uint_t N, double *wk, uint_t nwk, uint_t K, double *Sc);

// CONFIDENCE_FACTOR
// Computes the confidence factors s.t. confidence interval is S*Cf
// Cf is a vector of length nwk
double confidence_factor(CONF_BOUND side, double p, double v);
void confidence_factor(double p, double *wk, uint_t nwk, uint_t K, double *Cf);
void confidence_factor(double p, double *v, uint_t nv , double *Cf);

// FREQ_VALS
// fills the frequency vector f, values in range of 0 Fs
void freq_vals(double Fs, uint_t N, double *f);

//loads fftw if required
void fftw_load();
//cleans fftw if required
void fftw_clean();
