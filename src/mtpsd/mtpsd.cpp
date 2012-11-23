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

#include "mtpsd.h"
#include <cstdlib>      // needed for atexit to unload fftw properly
#include <cstring>      // needed for memcpy

/////////////////////////////////////////////////////////////////////
//  Template Instantiations
////////////////////////////////////////////////////////////////////

template  void eigenspecta<double>(double *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk, double *Sk);
template  void eigencoeffs<double>(double *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk);

template  void eigenspecta<fftw_complex>(fftw_complex *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk, double *Sk);
template  void eigencoeffs<fftw_complex>(fftw_complex *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk);

template  class mtpsd<double>; 
template  mtpsd<double>::mtpsd(double *data, uint_t n, double nW);
template  mtpsd<double>::mtpsd(double *data, mtpsd_workspace work);
template  mtpsd<double>::mtpsd(double *data, const double *tapers, const double *lambda, uint_t n, uint_t K);
template  mtpsd<double>::mtpsd(double *data, const double *tapers, const double *lambda, mtpsd_workspace work);
template  mtpsd<double>::~mtpsd();
template  void mtpsd<double>::compute();
template  double mtpsd<double>::operator()(uint_t ii);
template  double mtpsd<double>::operator()(uint_t K, uint_t ii);
template  fftw_complex mtpsd<double>::eig_coeff(uint_t K, uint_t ii);
template  double mtpsd<double>::freq(uint_t ii);
template  double mtpsd<double>::dof(uint_t ii);
template  double mtpsd<double>::wt(uint_t kk, uint_t ii);
template  double mtpsd<double>::conf_int(uint_t ii, CONF_BOUND side, double p);
template  double mtpsd<double>::conf_factor(uint_t ii, CONF_BOUND side, double p);
template  double mtpsd<double>::lambda(uint_t k);
template  double mtpsd<double>::taper(uint_t k, uint_t i);
template  double mtpsd<double>::F_stat(uint_t ii);
template  double mtpsd<double>::F_thresh(uint_t ii, double p);
template  bool mtpsd<double>::F_test(uint_t ii, double p);
template  uint_t mtpsd<double>::length();
template  uint_t mtpsd<double>::size(int dim);
template  const double *mtpsd<double>::pS();
template  mtpsd_workspace mtpsd<double>::getinfo();

template  class mtpsd<fftw_complex>; 
template  mtpsd<fftw_complex>::mtpsd(fftw_complex *data, uint_t n, double nW);
template  mtpsd<fftw_complex>::mtpsd(fftw_complex *data, mtpsd_workspace work);
template  mtpsd<fftw_complex>::mtpsd(fftw_complex *data, const double *tapers, const double *lambda, uint_t n, uint_t K);
template  mtpsd<fftw_complex>::mtpsd(fftw_complex *data, const double *tapers, const double *lambda, mtpsd_workspace work);
template  mtpsd<fftw_complex>::~mtpsd();
template  void mtpsd<fftw_complex>::compute();
template  double mtpsd<fftw_complex>::operator()(uint_t ii);
template  double mtpsd<fftw_complex>::operator()(uint_t K, uint_t ii);
template  fftw_complex mtpsd<fftw_complex>::eig_coeff(uint_t K, uint_t ii);
template  double mtpsd<fftw_complex>::freq(uint_t ii);
template  double mtpsd<fftw_complex>::dof(uint_t ii);
template  double mtpsd<fftw_complex>::wt(uint_t kk, uint_t ii);
template  double mtpsd<fftw_complex>::conf_int(uint_t ii, CONF_BOUND side, double p);
template  double mtpsd<fftw_complex>::conf_factor(uint_t ii, CONF_BOUND side, double p);
template  double mtpsd<fftw_complex>::lambda(uint_t k);
template  double mtpsd<fftw_complex>::taper(uint_t k, uint_t i);
template  double mtpsd<fftw_complex>::F_stat(uint_t ii);
template  double mtpsd<fftw_complex>::F_thresh(uint_t ii, double p);
template  bool mtpsd<fftw_complex>::F_test(uint_t ii, double p);
template  uint_t mtpsd<fftw_complex>::length();
template  uint_t mtpsd<fftw_complex>::size(int dim);
template  const double *mtpsd<fftw_complex>::pS();
template  mtpsd_workspace mtpsd<fftw_complex>::getinfo();

/////////////////////////////////////////////////////////////////////
//  MTPSD OBJECT
/////////////////////////////////////////////////////////////////////

void fix_workspace(mtpsd_workspace &work){
    
    if (work.n<0) work.n=0;
    
    if (work.nW<1 || work.nW > work.n/2.0) 
        work.nW=1;
    
    if ( work.K > work.n || work.K < 2)
        work.K=floor(2*work.nW)-1;
    work.K=tmath::max<uint_t>(work.K,2);       //require at least two sequences
    
    if (work.Fs<0) work.Fs=1;    
    
    work.N=tmath::max<uint_t>(work.N, work.n);
    
    if (work.weight_method==ADAPT)
        work.nwk=work.N;
    else
        work.nwk=1;

}

template <class T>
void mtpsd<T>::create_arrays(){
 
    //double check whether real data or complex data is entered
    if (sizeof(T)>sizeof(double)){
        this->info.dtype=COMPLEX_DATA;
    }
    else {
        this->info.dtype=REAL_DATA;
    }
    
    this->x = new T[this->info.n];
    this->S=new double[this->info.N];
    this->wk=new double[this->info.nwk * this->info.K];
    this->h=new double[this->info.n * this->info.K];
    this->l=new double[this->info.K];
    
    this->Jk=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->info.N * this->info.K);
    
}

template <class T>
void mtpsd<T>::init(T *data, uint_t n, uint_t K){
    
    this->info.weight_method=ADAPT;
    this->info.n=n;
    this->info.N=n;
    this->info.nwk=n;
    this->info.K=K;
    this->info.Fs=n;
    fix_workspace(this->info);
    
    this->create_arrays();
    memcpy(this->x, data, sizeof(T)*this->info.n);
    this->tapers_ready=false;
    
}

template <class T>
void mtpsd<T>::init(T *data, mtpsd_workspace work){
 
    this->info=work;
    fix_workspace(this->info);
    
    this->create_arrays();
    memcpy(this->x, data, sizeof(T)*this->info.n);
    this->tapers_ready=false;
    
}

template <class T>
mtpsd<T>::mtpsd(T *data, uint_t n, double nW){
    
    this->init(data,n, floor(2*nW)-1);
    this->info.nW=nW;
    
}

template <class T>
mtpsd<T>::mtpsd(T *data, mtpsd_workspace work){
    
    this->init(data,work);
    
}

template <class T>
mtpsd<T>::mtpsd(T *data, const double *tapers, const double *lambda, uint_t n, uint_t K){
    
    this->init(data,n,K);
    this->info.nW=-1;   //unknown
    
    //copy tapers
    memcpy(this->l, lambda, sizeof(double)*K);
    memcpy(this->h, tapers, sizeof(double)*this->info.n*K);
    
    if (K < 2){
        for (uint_t ii=K; ii<2; ii++){
            for (uint_t jj=0; jj<this->info.n; jj++){
                this->h[ii*this->info.n+jj]=0;  //zero out any other tapers
            }
            this->l[ii]=0;      //zero out eigenvalue
        }
    }
    
    this->tapers_ready=true;
    
}

template <class T>
mtpsd<T>::mtpsd(T *data, const double *tapers, const double *lambda, mtpsd_workspace work){
    
    this->init(data,work);
    
    memcpy(this->l, lambda, sizeof(double)*this->info.K);
    memcpy(this->h, tapers, sizeof(double)*this->info.n * this->info.K);
    this->tapers_ready=true;
    
}

template <class T>
mtpsd<T>::~mtpsd(){
    
    //cleanup
    fftw_free(this->Jk);
    delete [] this->x;    
    delete [] this->S;
    delete [] this->wk;
    delete [] this->h;
    delete [] this->l;
    
}

template <class T>
double mtpsd<T>::operator()(uint_t ii){
    return this->S[ii];
}

template <class T>
double mtpsd<T>::operator()(uint_t K, uint_t ii){
    double Sk=cabs( this->Jk[ K * this->info.N + ii ]);
    Sk=Sk*Sk;
    return Sk;
}

template <class T>
double mtpsd<T>::freq(uint_t ii){
    return (ii%this->info.N)*this->info.Fs/this->info.N;
}

template <class T>
fftw_complex mtpsd<T>::eig_coeff(uint_t K, uint_t ii){
    return Jk[ K* this->info.N + ii ];
}

template <class T>
void mtpsd<T>::compute(){
    
    if (this->tapers_ready==false){
        
        dpss_workspace taper_info;
        taper_info.n=this->info.n;
        taper_info.nW=this->info.nW;
        taper_info.seql=0;
        taper_info.sequ=this->info.K-1;
        taper_info.K=this->info.K;
        taper_info.interp_method=NONE;
        taper_info.interp_base=this->info.n;
        taper_info.energy=true;
        fix_workspace(taper_info);
        
        dpss mytapers(taper_info);
        
        try{ mytapers.compute(); }
        catch(...){  throw; }       //re-throws error
        
        memcpy(this->h, mytapers.ph(), sizeof(double)*taper_info.K*taper_info.n);
        memcpy(this->l, mytapers.pl(), sizeof(double)*taper_info.K);
        
        this->tapers_ready=true;
        
    }
    
    double varx;
    if (this->info.remove_mean == true)
        varx = tmath::var<T>(this->x,this->info.n);
    else
        varx = tmath::mom2<T>(this->x,this->info.n);
    
    //compute individual eigenspectra
    eigencoeffs<T>(this->x,this->info.n, this->h, this->l, this->info.K, this->info.remove_mean, this->info.N, this->Jk);
    
    //compute weights and final spectrum estimate
    if (this->info.weight_method==ADAPT){ 
        
        //construct adaptive weights and spectrum estimate
        adapt_wk(varx, this->l, this->Jk, this->info.K, this->info.N, ADAPT_TOL, this->wk, this->S);
    }
    else {
        
        //calculate fixed weights
        if (this->info.weight_method==EQUAL){
            for ( uint_t ii=0; ii<this->info.K; ii++)
                this->wk[ii]=1.0/this->info.K;
        }
        else{
            double lsum=tmath::sum<double>(this->l, this->info.K);
            for ( uint_t ii=0; ii<this->info.K; ii++)
                this->wk[ii]=this->l[ii]/lsum;            
        }
        
        //form final spectrum estimate
        combine_eig_coeffs( this->Jk, this->info.K, this->info.N, this->wk, this->info.nwk, this->S);
        
    }
    
    //scale for non-unit Fs
    tmath::scale<double>(this->S, 1.0/this->info.Fs,this->info.N, this->S);
    tmath::scale<fftw_complex>(this->Jk, 1.0/sqrt(this->info.Fs), this->info.N*this->info.K, this->Jk);
    
}


template <class T>
double mtpsd<T>::wt(uint_t kk,uint_t ii){
    
    return this->wk[ii%this->info.nwk+kk*this->info.nwk];
    
}

template <class T>
uint_t mtpsd<T>::length(){
    
    return this->info.N;
    
}

template <class T>
uint_t mtpsd<T>::size(int dim){

    if ( dim < 1 )
        return this->info.N;
    return this->info.K;
    
}

template <class T>
const double *mtpsd<T>::pS(){
        return this->S;
}

template <class T>
mtpsd_workspace mtpsd<T>::getinfo(){
        return this->info;
}

template <class T>
double mtpsd<T>::dof(uint_t ii){
 
   return degrees_of_freedom(ii,this->wk, this->info.nwk,this->info.K);
    
}

template <class T>
double mtpsd<T>::conf_int(uint_t ii, CONF_BOUND side, double p){

    return S[ii]*this->conf_factor(ii,side,p);
    
}

template <class T>
double mtpsd<T>::conf_factor(uint_t ii, CONF_BOUND side, double p){
    
    return confidence_factor(side,p,this->dof(ii));
    
}

template <class T>
double mtpsd<T>::F_stat(uint_t ii){
    
    return F_statistic(ii, this->Jk, this->info.N, this->info.K, this->wk, this->info.nwk, this->h, this->info.n);
    
}

template <class T>
double mtpsd<T>::F_thresh(uint_t ii, double p){
    
    double v= this->dof(ii)-2;  //have to subtract 2 for estimating eigencoefficient
    return F_threshold(p,v);
    
}

template <class T>
bool mtpsd<T>::F_test(uint_t ii, double p){
    
    double Fu=this->F_thresh(ii,p);
    double F=F_stat(ii);
    
    bool fail=true;
    
    if ( (Fu < 0) || (F < Fu) )
        fail=false;
    
    return fail;
    
}

template <class T>
double mtpsd<T>::lambda(uint_t k){
    
    if ( ( k < this->info.K) && (k >= 0) )
        return this->l[k];
    return 0;
    
}

template <class T>
double mtpsd<T>::taper(uint_t k, uint_t i){
    
    if ( ( k< this->info.K) && (k >= 0) && (i < this->info.n) && (i >= 0) )
        return this->h[ k*(this->info.n) + i];
    return 0;
    
}

/////////////////////////////////////////////////////////////////////
//  Code
/////////////////////////////////////////////////////////////////////

template <class T>
void eigencoeffs(T *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk){
        
    T m=tmath::mean<T>(x,n);       //computes standard mean
    
    double tmp;             //temporary variable
    T wm;                   //weighted mean
    const int nsize=N;
    
    // Remove weighted averages, reducing bias from non-centred data 
    //(forces eigenspectra to have zero DC component)
    for ( uint_t ii=0; ii<K; ii++){
        
        if (remove_mean==true){
            tmp=tmath::sum<double>(&tapers[ii*n],n);
        
            //tmp should be near zero for odd tapers (but due to round-off may not be exactly zero)
            if (  tmath::abs<double>(tmp) > ZERO_TOL ){
                wm=tmath::dot_mult<T,double>(x,&tapers[ii*n],n);
                wm=wm/tmp;
            }
            else {
                // for odd DPSS sequences, the weighted average is zero
                // However, we remove the regular mean for good measure
                wm=m;            
            } 
        }
        else {
            wm=0;
        }
 
        for (uint_t jj=0; jj<n; jj++){
            Jk[ii*N+jj]=x[jj]-wm;       //prepares for in-place FFT
        }        
    }
    
    //Window the data
    for ( uint_t ii=0; ii<K; ii++){
        tmath::pw_mult<fftw_complex,double>(&Jk[ii*N], &tapers[ii*n], n, &Jk[ii*N] );
        //zero-pad rest of array
        for (uint_t jj=n; jj<N; jj++){
            Jk[ii*N+jj]=0;
        }
        
    }
    
    // COMPUTE EIGENSPECTRA
    //computes all K ffts in one step
    fftw_load(); //potentially load fftw3
    fftw_plan p = fftw_plan_many_dft(1, &nsize, K , Jk, NULL, 1, N, Jk, NULL, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
 
}

// fills in the eigenspectum
template <class T>
void eigenspecta(T *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jk, double *Sk){
    
    eigencoeffs<T>(x,n,tapers,lambda,K,remove_mean,N,Jk);
    
    //compute eigenspectra
    double magJk;
    for (uint_t ii=0; ii<N*K; ii++){
        magJk=cabs(Jk[ii]);
        Sk[ii]= magJk*magJk;
    }
    
}

// Uses adaptive weighting to compute the weights for the eigenspectra, wk, and 
// fills the final spectrum estimate S.  The adaptive procedure iterates until 
// || S_new - S_old ||/N < tol
void adapt_wk( double varx, const double *lambda, const double *Sk, uint_t K, uint_t N, double tol, double *wk, double *S_init, double *S){
    
    uint_t nwk=N;        
    
    //copy initialization into estimate if required
    if (S!=S_init)  //this should be checking pointer values
        memcpy(S,S_init, N*sizeof(double));
    
    double diff = 2*tol;                        //forces while loop to start
    double b,b2, wdenom;                        //temporary variables
    while (diff > tol){
        
        for (uint_t ii=0; ii<N; ii++){
            
            //compute K weights corresponding to frequency ii
            wdenom=0;
            for ( uint_t kk=0; kk < K; kk++){
                b = S[ii] / (S[ii]*lambda[kk] + varx*(1-lambda[kk]) ); 
                b2 = b*b;
                wk[ii+kk*N]=b2*lambda[kk];
                wdenom+=wk[ii+kk*N];        //denominator for these K weights (to normalize)
            }
            
            // if all weights are zero (occurs of Sk[ii]=0 for all K)
            // then restore equal weighting, otherwise normalize
            if (wdenom==0){
                for ( uint_t jj=0; jj < K; jj++)
                    wk[ii+jj*N]=1/K;
            }
            else{
                for ( uint_t jj=0; jj < K; jj++)
                    wk[ii+jj*N]=wk[ii+jj*N]/wdenom;   //normalize weights
            }           
        }
        
        //computes new estimate using the adapted weights, 
        // AND calculates the normalized difference
        combine_eig_spec( Sk, K, N, wk, nwk, S, diff);                
    }
    
}

void adapt_wk( double varx, const double *lambda, const fftw_complex *Jk, uint_t K, uint_t N, double tol, double *wk, double *S_init, double *S){
    
    uint_t nwk=N;        
    
    //copy initialization into estimate if required
    if (S!=S_init)  //this should be checking pointer values
        memcpy(S,S_init, N*sizeof(double));
    
    double diff = 2*tol;                        //forces while loop to start
    double b,b2, wdenom;                        //temporary variables
    
    while (diff > tol){
        
        for (uint_t ii=0; ii<N; ii++){
            
            //compute K weights corresponding to frequency ii
            wdenom=0;
            for ( uint_t kk=0; kk < K; kk++){
                b = S[ii] / (S[ii]*lambda[kk] + varx*(1-lambda[kk]) ); 
                b2 = b*b;
                wk[ii+kk*N]=b2*lambda[kk];
                wdenom+=wk[ii+kk*N];        //denominator for these K weights (to normalize)
            }
            
            // if all weights are zero (occurs of Sk[ii]=0 for all K)
            // then restore equal weighting, otherwise normalize
            if (wdenom==0){
                for ( uint_t jj=0; jj < K; jj++)
                    wk[ii+jj*N]=1/K;
            }
            else{
                for ( uint_t jj=0; jj < K; jj++)
                    wk[ii+jj*N]=wk[ii+jj*N]/wdenom;   //normalize weights
            }           
        }
        
        //computes new estimate using the adapted weights, 
        // AND calculates the normalized difference
        combine_eig_coeffs( Jk, K, N, wk, nwk, S, diff);
        
    }
    
}

void adapt_wk( double varx, const double *lambda, const double *Sk, uint_t K, uint_t N, double tol, double *wk, double *S){
    
    //first, set up weights to find average of first two Sk
    for ( uint_t ii=0; ii<K; ii++)
        wk[ii]=0;
    wk[0]=0.5;
    wk[1]=0.5;
    combine_eig_spec( Sk, 2, N, wk, 1, S);      //S = 1/2 ( S_0 + S_1 );
    adapt_wk(varx,lambda,Sk,K,N,tol,wk,S,S);      // uses current S as initialization
    
}

void adapt_wk( double varx, const double *lambda, const fftw_complex *Jk, uint_t K, uint_t N, double tol, double *wk, double *S){
    
    //first, set up weights to find average of first two Sk
    for ( uint_t ii=0; ii<K; ii++)
        wk[ii]=0;
    wk[0]=0.5;
    wk[1]=0.5;
    combine_eig_coeffs( Jk, 2, N, wk, 1, S);      //S = 1/2 ( S_0 + S_1 );
    
    adapt_wk(varx,lambda,Jk,K,N,tol,wk,S,S);      // uses current S as initialization
    
}

void combine_eig_spec( const double *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S){
    
    uint_t iwt; //index of weight to use
    
    for (uint_t ii=0; ii<N; ii++){
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk[jj*N+ii];   //weighted average
        }
    }
    
}

void combine_eig_coeffs( const fftw_complex *Jk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S){
    
    uint_t iwt; //index of weight to use
    double Sk=0;
    for (uint_t ii=0; ii<N; ii++){
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            Sk=cabs(Jk[jj*N+ii]);
            Sk=Sk*Sk;
            
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk;   //weighted average
        }
    }
    
}

// NOTE: S is both an input AND output.
// Combines eigenspectra, fills S, but also computes the normalized absolute difference from what
// was previously in S.  This difference = || S_new - S_old ||/N is returned in diff.
void combine_eig_spec( const double *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S, double &diff){
    
    uint_t iwt; //index of weight to use
    
    diff=0;
    double Sl, absdiff;
    for (uint_t ii=0; ii<N; ii++){
        Sl=S[ii];                               //stores previous value
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk[jj*N+ii];   //weighted average
        }
        absdiff=tmath::abs<double>(S[ii]-Sl);          //absolute difference
        diff=diff+absdiff;
    }
    
    diff=sqrt(diff)/N;
    
}

void combine_eig_coeffs( const fftw_complex *Jk, uint_t K, uint_t N, const double *wk, uint_t nwk, double *S, double &diff){
    
    uint_t iwt; //index of weight to use
    double Sk;
    
    diff=0;
    double Sl, absdiff;
    for (uint_t ii=0; ii<N; ii++){
        Sl=S[ii];                               //stores previous value
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            Sk=cabs(Jk[jj*N+ii]);
            Sk=Sk*Sk;
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk;   //weighted average
        }
        absdiff=tmath::abs<double>(S[ii]-Sl);          //absolute difference
        diff=diff+absdiff;
    }
    
    diff=sqrt(diff)/N;
    
}


//should be bounded below by 2 and above by 2K (if wk are normalized correctly)
double degrees_of_freedom(uint_t ii, double *wk, uint_t nwk, uint_t K){
 
    double v=0;   // approximate degrees of freedom
    for ( uint_t kk=0; kk < K; kk++){
        v += wk[kk*nwk + ii%nwk]*wk[kk*nwk + ii%nwk];
    }
    v=2.0/v;   //should be bounded below by 2 and above by 2K (if wk are normalized correctly)
    return v;
    
}

// v is the equivalent degrees of freedom
double confidence_factor(CONF_BOUND side, double p, double v){
  
    if ( ( p <= 0 ) || ( p >= 100 ) )
        p=0;
    else if ( p >= 1 )
        p=p/100;
    
    double p_in;
    if (side==LOWER)
        p_in=(1.0-p)/2.0;
    else
        p_in=(1.0+p)/2.0;
    
    double Qc=stats::chisquared_inv(p_in, v);      //quantile
    
    return v/Qc;
    
}

double F_statistic(uint_t ii, fftw_complex *Jk, uint_t N, uint_t K, double *wk, uint_t nwk, double *h, uint_t n ){
 
    double ak, aH02=0, a2H02=0, sige2, Jkerr, wssq, tmp;
    fftw_complex C=0, Jke;     // estimate of fourier coefficient at f (assuming line frequency)
    
    double *H0=new double[K];
    
    for ( uint_t jj=0; jj< K; jj++){
        ak = sqrt( wk[jj*nwk + ii%nwk] );       // preserves relative weights in estimate of S
        H0[jj] = tmath::sum<double>( &h[jj*n], n ) ;
        if (tmath::abs<double>(H0[jj])<ZERO_TOL) H0[jj]=0;     // zeroes out small values
        tmp=H0[jj]*H0[jj]*ak;
        aH02 += tmp;
        a2H02+= tmp*ak;
        C += Jk[ii+jj*N]*ak*H0[jj];
    }
    
    C=C/aH02;   // estimated fourier coefficient
    sige2=0;
    wssq=0;     // sum of squares of wk (used in estamating equiv. DOF)
    uint_t iwk;
    for ( uint_t jj=0; jj< K; jj++){
        iwk = jj*nwk + ii%nwk;
        //H0 = sum<double>( &h[jj*n], n ) ;      //re-calculates in order to prevent creating an array
        Jke = C*H0[jj];
        Jkerr = cabs(Jk[ii+jj*N]-Jke);
        sige2 += wk[iwk]*Jkerr*Jkerr;   //weighted error
        wssq += wk[iwk]*wk[iwk];
    }
    
    //clean up before return
    delete [] H0;
    
    // The two independent chi^2  distributions for the modified F-test are:
    //
    //   [ ( 2|C|^2 )*( sum( ak*H(0)^2,   k=0..K-1) )^2 ] /      ~ chi^2_2
    //   [  sige^2    ( sum( ak^2*H(0)^2, k=0..K-1) )   ]
    //
    //   [  2 sigehat^2 ] / [  sige^2 wssq ]      ~ chi^2_( 2/wssq - 2 )
    //
    // NOTE: This F-test will fail at any frequency where the estimate
    //       has near two degrees of freedom (most of the weight is on
    //       a single eigenspectrum).  This MAY occur with adaptive weighting
    //       if the estimate at f is dominated by bias.  In such cases, 
    //       either the F-test should be performed with equal (or near-equal),
    //       or it should be concluded there is no significant frequency
    //       content at f (since it is drowned out by bias anyway)
    
    return ( cabs(C)*cabs(C) ) * (1-wssq) * aH02*aH02 / sige2 / a2H02;
    
}

double F_threshold(double p, double v){
    
    if ( ( p <= 0 ) || ( p >= 100 ) )
        p=0;
    else if ( p >= 1 )
        p=p/100;
    
    double a=1-p;
    double b=pow(a,2.0/v);
    double Fu=INFINITY;   //defaults to division by zero
    
    if ( b > ZERO_TOL )
        Fu=v*(1.0-b)/2.0/b;
    
    return Fu;
    
}

// Calculates the equivalent degrees of freedom of the estimate (assuming large n)
// such that S_estimate ~ S_true chi^2_v / v,
void degrees_of_freedom( double *wk, uint_t nwk, uint_t K, double *v){
    
    //fills v with 2/sum_k(wk^2)  [ = 2K for unity weights ]
    for (uint_t ii=0; ii<nwk; ii++){
        v[ii]=degrees_of_freedom( ii, wk, nwk, K);   
    }
    
}

void confidence_interval( double p, double *S, uint_t N, double *v, uint_t nv, double *Sc){
 
    uint_t pp_start=N-nv;       //temporarily stores quantile information in Sc (ensuring not to be overwritten)
    
    stats::chisquared_inv((1.0-p)/2.0, v, nv, &Sc[pp_start]);   //lower quantile
    stats::chisquared_inv((1.0+p)/2.0, v, nv, &Sc[pp_start+N]); //upper quantile
    
    uint_t ipp;                 //used to calculate index of quantile
    
    for (uint_t ii=0; ii<N; ii++){
        ipp=ii%nv;
        Sc[ii]   = v[ipp]*S[ii]/Sc[pp_start+ipp];   //lower bound
        Sc[ii+N] = v[ipp]*S[ii]/Sc[pp_start+N+ipp]; //upper bound
    }
}

void confidence_interval( double p, double *S, uint_t N, double *wk, uint_t nwk, uint_t K, double *Sc){
    
    uint_t pp_start=N-nwk;       //temporarily stores quantile information in Sc (ensuring not to be overwritten)
    
    double v;   // approximate degrees of freedom
    for (uint_t ii=0; ii<nwk; ii++){
        v=degrees_of_freedom(ii,wk,nwk,K);        
        Sc[pp_start+ii]=confidence_factor(LOWER,p,v);
        Sc[pp_start+N+ii]=confidence_factor(UPPER,p,v);
        
    }
    
    uint_t ipp;                 //used to calculate index of quantile
    for (uint_t ii=0; ii<N; ii++){
        ipp=ii%nwk;
        Sc[ii]   = Sc[ipp+pp_start]*S[ii];          //lower bound
        Sc[ii+N] = Sc[ipp+pp_start+N]*S[ii];        //upper bound
    }
}

// computes the confidence factors
void confidence_factor(double p, double *wk, uint_t nwk, uint_t K, double *Cf){
    
    double v;   // approximate degrees of freedom
    for (uint_t ii=0; ii<nwk; ii++){
        v=degrees_of_freedom(ii,wk,nwk,K);        
        Cf[ii]=confidence_factor(LOWER,p,v);
        Cf[ii+nwk]=confidence_factor(UPPER,p,v);
        
    }
    
}

void confidence_factor(double p, double *v, uint_t nv , double *Cf){
    
    for (uint_t ii=0; ii<nv; ii++){
        Cf[ii]=confidence_factor(LOWER,p,v[ii]);
        Cf[ii+nv]=confidence_factor(UPPER,p,v[ii]);
    }
    
}

// computes F-statistic
void F_statistic(fftw_complex *Jk, uint_t N, uint_t K, double *wk, uint_t nwk, double *h, uint_t n, double *F){
    
    for (uint_t ii=0; ii<N; ii++){
        F[ii] = F_statistic(ii,Jk, N, K, wk, nwk, h, n);
    }
    
}

void F_threshold(double p, double *wk, uint_t nwk, uint_t K, double *Fu){
    
    double v;
    for (uint_t ii=0; ii<nwk; ii++){
        v=degrees_of_freedom(ii,wk,nwk,K);
        Fu[ii] = F_threshold(p,v);
    }
    
}

void F_threshold(double p, double *v, uint_t nv, double *Fu){
    
    for (uint_t ii=0; ii<nv; ii++){
        Fu[ii] = F_threshold(p,v[ii]);
    }
    
}

// fills the frequency vector f, values in range of 0 Fs
void freq_vals(double Fs, uint_t N, double *f){
    
    double df=1.0/N*Fs;
    for (uint_t ii=0; ii<N; ii++){
        f[ii]=ii*df;
    }
    
}

void fftw_load(){
    
    static bool fftw_loaded = false;
    if (fftw_loaded==false){
        fftw_loaded=true;
        atexit(fftw_clean); // if fftw is not already loaded, then set cleanup function on exit
    }
    
}

void fftw_clean(){
    //cleanup routines
    fftw_cleanup();
}

