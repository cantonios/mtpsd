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

/////////////////////////////////////////////////////////////////////////////////////////////
//  mtpsd.cpp contains all functions required for calculating the multitaper
//            spectrum estimate for a one-dimensional time-series (real or complex)
//            
//            An interface class is provided ( class mtpsd{} ) to simplify use.
//            The individual routines are also supplied that work directly on 
//            arrays if you wish to have more control over memory use.
//
/////////////////////////////////////////////////////////////////////////////////////////////

#define MTPSD_VERSION 1.0

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

enum WEIGHT_METHOD {ADAPT, EIGEN, EQUAL};
enum DATA_TYPE {REAL_DATA, COMPLEX_DATA};
enum CONF_BOUND {LOWER, UPPER};

// Constants involved in computing psd
struct mtpsd_workspace{    
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
    mtpsd_workspace(): n(0), nW(2), weight_method(ADAPT), N(0), K(3), Fs(1), remove_mean(true), dtype(REAL_DATA), nwk(0){ }
    mtpsd_workspace(uint_t n_, double nW_):n(n_),nW(nW_), weight_method(ADAPT), N(n_),K(floor(2*nW_)-1), Fs(1), remove_mean(true), dtype(REAL_DATA), nwk(n_){ }
};

void fix_workspace(mtpsd_workspace &work);      //fixes possibly invalid entries

template <class T>
class mtpsd{
    
public:    
    mtpsd(T *data, uint_t n, double nW);     //uses dpss{} class to compute tapers
    mtpsd(T *data, mtpsd_workspace work);
    mtpsd(T *data, const double *tapers, const double *lambda, uint_t n, uint_t K);      //allows k custom tapers of length n
    mtpsd(T *data, const double *tapers, const double *lambda, mtpsd_workspace work);
          
    ~mtpsd();
    
    void compute();                          // Computes Spectrum with current information
    
    double operator()(uint_t ii);            // returns the ith value of the spectrum estimate
    double operator()(uint_t k, uint_t ii);  // returns the ith value of eigenspectrum k
    fftw_complex eig_coeff(uint_t k, uint_t ii);   //returns the ith value of eigencoefficient k (Jk) 
    
    double freq(uint_t ii);                 // returns ith frequency
    
    double dof(uint_t ii);                  // returns equivalent degrees of freedom of estimate at frequency ii
    double wt(uint_t kk, uint_t ii);        // returns w_kk(ii), the kkth weight at frequency ii
    double conf_int(uint_t ii, CONF_BOUND side, double p);     //returns either lower or upper bound of px100% confidence interval at frequency ii
    double conf_factor(uint_t ii, CONF_BOUND side, double p);  //returns either lower or upper confidence factor (conf_int=S[ii]*conf_factor)
    
    double lambda(uint_t kk);                //returns the kth energy concentration, seql<= k <= sequ
    double taper(uint_t kk, uint_t ii);      //returns the iith value of sequence seql<= k <= sequ
    
    double F_stat(uint_t ii);               //returns F-statistic of frequency ii
    double F_thresh(uint_t ii, double p);   //returns the F-statistic threshold at frequency ii with px100% acceptance probability
    bool F_test(uint_t ii, double p);       //returns result of F_test, true if frequency ii is significant
    
    uint_t length();                        // returns N, length of S
    uint_t size(int dim);                   //if dim < 1 returns length of S, otherwise returns K
    const double *pS();                     //returns pointer to Spectrum S
    
    mtpsd_workspace getinfo();              //returns information about this mtpsd class
      
protected:
    
    void init(T *data, uint_t n, uint_t K);     //initializes variables
    void init(T *data, mtpsd_workspace work);
    void create_arrays();                       //allocates memory
    
    bool tapers_ready;      // true if tapers are already computed or supplied, false otherwise
    
    T *x;                   // Original Data
    double *S;              // Final spectrum
    fftw_complex *Jk;       // Eigencoefficients, |Jk|^2 = Sk
    double *wk;             // weights (unity, energy concentrations, or adaptive);
    double *h;              // data tapers
    double *l;              // taper energies
    mtpsd_workspace info;   // holds information about the PSD
    
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
void eigenspecta(T *x, uint_t n, const double *h, const double *l, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk, double *Sk);
template <class T>
void eigencoeffs(T *x, uint_t n, const double *h, const double *l, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk);

// COMBINE_EIGSPEC
// combine eigenspectra using weights wk, with size (K x nwk)
// nwk is the # of unique weights per eigenspectra, is usually nwk = 1 (same weight across all freq) or N (different across all f)
// Weights must be normalized so they sum to one at each frequency
void combine_eig_spec( const double *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S);
void combine_eig_coeffs( const fftw_complex *Jk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S);

// Combines eigenspectra, fills S, but also computes the normalized absolute difference from
// the previous.  This difference = || S_new - S_old ||/N is returned in diff.
// NOTE: S is both an input AND output in this procedure.
void combine_eig_spec( const double *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S, double &diff);
void combine_eig_coeffs( const fftw_complex *Jk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S, double &diff);

// Uses adaptive weighting to compute the weights for the eigenspectra, wk, and 
// fills the final spectrum estimate S.  The adaptive procedure iterates until 
// || S_new - S_old ||/N < tol
void adapt_wk( double varx, const double *l, const double *Sk, uint_t K, uint_t N, double tol, double *wk, double *S_init, double *S);
void adapt_wk( double varx, const double *l, const fftw_complex *Jk, uint_t K, uint_t N, double tol, double *wk, double *S_init, double *S);
void adapt_wk( double varx, const double *l, const double *Sk, uint_t K, uint_t N, double tol, double *wk, double *S);
void adapt_wk( double varx, const double *l, const fftw_complex *Jk, uint_t K, uint_t N, double tol, double *wk, double *S);

// DEGREES_OF_FREEDOM
// Calculates the equivalent degrees of freedom of the estimate (assuming large n)
// such that S_estimate ~ S_true chi^2_v / v,
double degrees_of_freedom(uint_t ii, double *wk, uint_t nwk, uint_t K);     // at frequency ii
void degrees_of_freedom( double *wk, uint_t nwk, uint_t K, double *v);      // fills vector v

// CONFIDENCE_INTERVAL
// Computes the px100% confidence interval of the spectrum estimate S
// v is the degrees of freedom (vector of length nv)  (nv = 1 for constant dof across specturm, nv=N for varying)
// Sc has length 2*N, and stores the lower then upper confidence bounds in row-major format.
void confidence_interval( double p, double *S, uint_t N, double *v, uint_t nv, double *Sc);
void confidence_interval( double p, double *S, uint_t N, double *wk, uint_t nwk, uint_t K, double *Sc);

// CONFIDENCE_FACTOR
// Computes the confidence factors s.t. confidence interval is S*Cf
// Cf is a vector of length nwk
double confidence_factor(CONF_BOUND side, double p, double v);
void confidence_factor(double p, double *wk, uint_t nwk, uint_t K, double *Cf);
void confidence_factor(double p, double *v, uint_t nv , double *Cf);

// F_STATISTIC
// Returns the F-test statistic of the spectrum estimate based on the eigencoefficients (Jk), 
// the tapers (h), and the weights for the final estimate (wk). 
// Fills F-statistic array F
// NOTE: This F-test is modified from Percival/Walden (1998) pp. 499 to include 
//       weights (either adaptive or energy concentrations).
double F_statistic(uint_t ii, fftw_complex *Jk, uint_t N, uint_t K, double *wk, uint_t nwk, double *h, uint_t n );  // at freq ii
void F_statistic(fftw_complex *Jk, uint_t N, uint_t K, double *wk, uint_t nwk, double *h, uint_t n, double *F);

// F_THRESHOLD
// Computes the upper threshold for the F-test with px100% acceptance probability
// Fills Fu, length nwk
double F_threshold(double p, double v);
void F_threshold(double p, double *wk, uint_t nwk, uint_t K, double *Fu);
void F_threshold(double p, double *v, uint_t nv, double *Fu);

// FREQ_VALS
// fills the frequency vector f, values in range of 0 Fs
void freq_vals(double Fs, uint_t N, double *f);

//loads fftw if required
void fftw_load();
//cleans fftw if required
void fftw_clean();
