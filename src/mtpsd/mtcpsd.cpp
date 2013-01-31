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

#include "mtcpsd.h"
#include <cstdlib>      // needed for atexit to unload fftw properly
#include <cstring>      // needed for memcpy

/////////////////////////////////////////////////////////////////////
//  Template Instantiations
////////////////////////////////////////////////////////////////////

template  void eigenspecta<double>(double *x, double *y, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx, fftw_complex *Jky, fftw_complex *Sk);
template  void eigencoeffs<double>(double *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx);

template  void eigenspecta<fftw_complex>(fftw_complex *x, fftw_complex *y, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx, fftw_complex *Jky, fftw_complex *Sk);
template  void eigencoeffs<fftw_complex>(fftw_complex *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx);

template  class mtcpsd<double>;
template  mtcpsd<double>::mtcpsd(double *x, double *y, uint_t n, double nW);
template  mtcpsd<double>::mtcpsd(double *x, double *y, mtcpsd_workspace work);
template  mtcpsd<double>::mtcpsd(double *x, double *y, const double *tapers, const double *lambda, uint_t n, uint_t K);
template  mtcpsd<double>::mtcpsd(double *x, double *y, const double *tapers, const double *lambda, mtcpsd_workspace work);
template  mtcpsd<double>::~mtcpsd();
template  void mtcpsd<double>::compute();
template  fftw_complex mtcpsd<double>::operator()(uint_t ii);
template  fftw_complex mtcpsd<double>::operator()(uint_t K, uint_t ii);
template  fftw_complex mtcpsd<double>::eig_coeff_x(uint_t K, uint_t ii);
template  fftw_complex mtcpsd<double>::eig_coeff_y(uint_t K, uint_t ii);
template  double mtcpsd<double>::freq(uint_t ii);
template  double mtcpsd<double>::dof(uint_t ii);
template  double mtcpsd<double>::wt(uint_t kk, uint_t ii);
template  double mtcpsd<double>::conf_int(uint_t ii, CONF_BOUND side, double p);
template  double mtcpsd<double>::conf_factor(uint_t ii, CONF_BOUND side, double p);
template  double mtcpsd<double>::lambda(uint_t k);
template  double mtcpsd<double>::taper(uint_t k, uint_t i);
template  uint_t mtcpsd<double>::length();
template  uint_t mtcpsd<double>::size(int dim);
template  const fftw_complex *mtcpsd<double>::pS();
template  mtcpsd_workspace mtcpsd<double>::getinfo();

template  class mtcpsd<fftw_complex>;
template  mtcpsd<fftw_complex>::mtcpsd(fftw_complex *x, fftw_complex *y, uint_t n, double nW);
template  mtcpsd<fftw_complex>::mtcpsd(fftw_complex *x, fftw_complex *y, mtcpsd_workspace work);
template  mtcpsd<fftw_complex>::mtcpsd(fftw_complex *x, fftw_complex *y, const double *tapers, const double *lambda, uint_t n, uint_t K);
template  mtcpsd<fftw_complex>::mtcpsd(fftw_complex *x, fftw_complex *y, const double *tapers, const double *lambda, mtcpsd_workspace work);
template  mtcpsd<fftw_complex>::~mtcpsd();
template  void mtcpsd<fftw_complex>::compute();
template  fftw_complex mtcpsd<fftw_complex>::operator()(uint_t ii);
template  fftw_complex mtcpsd<fftw_complex>::operator()(uint_t K, uint_t ii);
template  fftw_complex mtcpsd<fftw_complex>::eig_coeff_x(uint_t K, uint_t ii);
template  fftw_complex mtcpsd<fftw_complex>::eig_coeff_y(uint_t K, uint_t ii);
template  double mtcpsd<fftw_complex>::freq(uint_t ii);
template  double mtcpsd<fftw_complex>::dof(uint_t ii);
template  double mtcpsd<fftw_complex>::wt(uint_t kk, uint_t ii);
template  double mtcpsd<fftw_complex>::conf_int(uint_t ii, CONF_BOUND side, double p);
template  double mtcpsd<fftw_complex>::conf_factor(uint_t ii, CONF_BOUND side, double p);
template  double mtcpsd<fftw_complex>::lambda(uint_t k);
template  double mtcpsd<fftw_complex>::taper(uint_t k, uint_t i);
template  uint_t mtcpsd<fftw_complex>::length();
template  uint_t mtcpsd<fftw_complex>::size(int dim);
template  const fftw_complex *mtcpsd<fftw_complex>::pS();
template  mtcpsd_workspace mtcpsd<fftw_complex>::getinfo();

/////////////////////////////////////////////////////////////////////
//  mtcpsd OBJECT
/////////////////////////////////////////////////////////////////////

void fix_workspace(mtcpsd_workspace &work){
    
    if (work.n<0) work.n=0;
    
    if (work.nW<1 || work.nW > work.n/2.0) 
        work.nW=1;
    
    if ( work.K > work.n || work.K < 2)
        work.K=floor(2*work.nW)-1;
    work.K=tmath::max<uint_t>(work.K,2);       //require at least two sequences
    
    if (work.Fs<0) work.Fs=1;    
    
    work.N=tmath::max<uint_t>(work.N, work.n);
    work.nwk=1;

}

template <class T>
void mtcpsd<T>::create_arrays(){
 
    //double check whether real data or complex data is entered
    if (sizeof(T)>sizeof(double)){
        this->info.dtype=COMPLEX_DATA;
    }
    else {
        this->info.dtype=REAL_DATA;
    }
    
    this->x = new T[this->info.n];
    this->y = new T[this->info.n];
    this->wk=new double[this->info.nwk * this->info.K];
    this->h=new double[this->info.n * this->info.K];
    this->l=new double[this->info.K];
    
    this->S=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*this->info.N);
    this->Jkx=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->info.N * this->info.K);
    this->Jky=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->info.N * this->info.K);
    
}

template <class T>
void mtcpsd<T>::init(T *x, T *y, uint_t n, uint_t K){
    
    this->info.weight_method=EIGEN;
    this->info.n=n;
    this->info.N=n;
    this->info.nwk=n;
    this->info.K=K;
    this->info.Fs=n;
    fix_workspace(this->info);
    
    this->create_arrays();
    memcpy(this->x, x, sizeof(T)*this->info.n);
    memcpy(this->y, y, sizeof(T)*this->info.n);
    this->tapers_ready=false;
    
}

template <class T>
void mtcpsd<T>::init(T *x, T *y,mtcpsd_workspace work){
 
    this->info=work;
    fix_workspace(this->info);
    
    this->create_arrays();
    memcpy(this->x, x, sizeof(T)*this->info.n);
    memcpy(this->y, y, sizeof(T)*this->info.n);
    this->tapers_ready=false;
    
}

template <class T>
mtcpsd<T>::mtcpsd(T *x, T* y, uint_t n, double nW){
    
    this->init(x, y, n, floor(2*nW)-1);
    this->info.nW=nW;
    
}

template <class T>
mtcpsd<T>::mtcpsd(T *x, T *y, mtcpsd_workspace work){
    
    this->init(x,y,work);
    
}

template <class T>
mtcpsd<T>::mtcpsd(T *x, T *y, const double *tapers, const double *lambda, uint_t n, uint_t K){
    
    this->init(x,y,n,K);
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
mtcpsd<T>::mtcpsd(T *x, T *y, const double *tapers, const double *lambda, mtcpsd_workspace work){
    
    this->init(x,y,work);
    
    memcpy(this->l, lambda, sizeof(double)*this->info.K);
    memcpy(this->h, tapers, sizeof(double)*this->info.n * this->info.K);
    this->tapers_ready=true;
    
}

template <class T>
mtcpsd<T>::~mtcpsd(){
    
    //cleanup
    fftw_free(this->Jkx);
    fftw_free(this->Jky);
    fftw_free(this->S);
    delete [] this->x;    
    delete [] this->wk;
    delete [] this->h;
    delete [] this->l;
    
}

template <class T>
fftw_complex mtcpsd<T>::operator()(uint_t ii){
    return this->S[ii];
}

template <class T>
fftw_complex mtcpsd<T>::operator()(uint_t K, uint_t ii){
    fftw_complex Sk= conj(this->Jkx[ K * this->info.N + ii ])*this->Jky[ K * this->info.N + ii ];
    return Sk;
}

template <class T>
double mtcpsd<T>::freq(uint_t ii){
    return (ii%this->info.N)*this->info.Fs/this->info.N;
}

template <class T>
fftw_complex mtcpsd<T>::eig_coeff_x(uint_t K, uint_t ii){
    return Jkx[ K* this->info.N + ii ];
}

template <class T>
fftw_complex mtcpsd<T>::eig_coeff_y(uint_t K, uint_t ii){
    return Jky[ K* this->info.N + ii ];
}

template <class T>
void mtcpsd<T>::compute(){
    
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
    
//    double varx;
//    if (this->info.remove_mean == true)
//        varx = tmath::var<T>(this->x,this->info.n);
//    else
//        varx = tmath::mom2<T>(this->x,this->info.n);
//
//    double vary;
//    if (this->info.remove_mean == true)
//		vary = tmath::var<T>(this->y,this->info.n);
//	else
//		vary = tmath::mom2<T>(this->y,this->info.n);

    //compute individual eigenspectra
    eigencoeffs<T>(this->x,this->info.n, this->h, this->l, this->info.K, this->info.remove_mean, this->info.N, this->Jkx);
    eigencoeffs<T>(this->y,this->info.n, this->h, this->l, this->info.K, this->info.remove_mean, this->info.N, this->Jky);
    
    //compute weights and final spectrum estimate
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
	combine_eig_coeffs( this->Jkx, this->Jky, this->info.K, this->info.N, this->wk, this->info.nwk, this->S);
    
    //scale for non-unit Fs
    tmath::scale<fftw_complex>(this->S,   1.0/this->info.Fs,this->info.N, this->S);
    tmath::scale<fftw_complex>(this->Jkx, 1.0/sqrt(this->info.Fs), this->info.N*this->info.K, this->Jkx);
    tmath::scale<fftw_complex>(this->Jky, 1.0/sqrt(this->info.Fs), this->info.N*this->info.K, this->Jky);
    
}


template <class T>
double mtcpsd<T>::wt(uint_t kk,uint_t ii){
    
    return this->wk[ii%this->info.nwk+kk*this->info.nwk];
    
}

template <class T>
uint_t mtcpsd<T>::length(){
    
    return this->info.N;
    
}

template <class T>
uint_t mtcpsd<T>::size(int dim){

    if ( dim < 1 )
        return this->info.N;
    return this->info.K;
    
}

template <class T>
const fftw_complex *mtcpsd<T>::pS(){
        return this->S;
}

template <class T>
mtcpsd_workspace mtcpsd<T>::getinfo(){
        return this->info;
}

template <class T>
double mtcpsd<T>::dof(uint_t ii){
 
   return degrees_of_freedom(ii,this->wk, this->info.nwk,this->info.K);
    
}

template <class T>
double mtcpsd<T>::conf_int(uint_t ii, CONF_BOUND side, double p){

    return cabs(S[ii])*this->conf_factor(ii,side,p);
    
}

template <class T>
double mtcpsd<T>::conf_factor(uint_t ii, CONF_BOUND side, double p){
    
    return confidence_factor(side,p,this->dof(ii));
    
}

template <class T>
double mtcpsd<T>::lambda(uint_t k){
    
    if ( ( k < this->info.K) && (k >= 0) )
        return this->l[k];
    return 0;
    
}

template <class T>
double mtcpsd<T>::taper(uint_t k, uint_t i){
    
    if ( ( k< this->info.K) && (k >= 0) && (i < this->info.n) && (i >= 0) )
        return this->h[ k*(this->info.n) + i];
    return 0;
    
}

/////////////////////////////////////////////////////////////////////
//  Code
/////////////////////////////////////////////////////////////////////

template <class T>
void eigencoeffs(T *x, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx){
        
    T mx=tmath::mean<T>(x,n);       //computes standard mean
    
    double tmp;             //temporary variable
    T wmx;                   //weighted mean
    const int nsize=N;
    
    // Remove weighted averages, reducing bias from non-centred data 
    //(forces eigenspectra to have zero DC component)
    for ( uint_t ii=0; ii<K; ii++){
        
        if (remove_mean==true){
            tmp=tmath::sum<double>(&tapers[ii*n],n);
        
            //tmp should be near zero for odd tapers (but due to round-off may not be exactly zero)
            if (  tmath::abs<double>(tmp) > ZERO_TOL ){
                wmx=tmath::dot_mult<T,double>(x,&tapers[ii*n],n);
                wmx=wmx/tmp;
            }
            else {
                // for odd DPSS sequences, the weighted average is zero
                // However, we remove the regular mean for good measure
                wmx=mx;
            } 
        }
        else {
            wmx=0;
        }
 
        for (uint_t jj=0; jj<n; jj++){
            Jkx[ii*N+jj]=x[jj]-wmx;       //prepares for in-place FFT
        }        
    }
    
    //Window the data
    for ( uint_t ii=0; ii<K; ii++){
        tmath::pw_mult<fftw_complex,double>(&Jkx[ii*N], &tapers[ii*n], n, &Jkx[ii*N] );
        //zero-pad rest of array
        for (uint_t jj=n; jj<N; jj++){
            Jkx[ii*N+jj]=0;
        }
    }
    
    // COMPUTE EIGENSPECTRA
    // computes all K ffts in one step
    fftw_load(); //potentially load fftw3
    fftw_plan px = fftw_plan_many_dft(1, &nsize, K , Jkx, NULL, 1, N, Jkx, NULL, 1, N, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(px);
    fftw_destroy_plan(px);
 
}

// fills in the eigenspectum
template <class T>
void eigenspecta(T *x, T *y, uint_t n, const double *tapers, const double *lambda, uint_t K, bool remove_mean, uint_t N, fftw_complex *Jkx, fftw_complex *Jky, fftw_complex *Sk){
    
    eigencoeffs<T>(x,n,tapers,lambda,K,remove_mean,N,Jkx);
    eigencoeffs<T>(y,n,tapers,lambda,K,remove_mean,N,Jky);
    
    //compute eigenspectra
    for (uint_t ii=0; ii<N*K; ii++){
        Sk[ii]= conj(Jkx[ii])*Jky[ii];
    }
    
}

void combine_eig_spec( const fftw_complex *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S){
    
    uint_t iwt; //index of weight to use
    
    for (uint_t ii=0; ii<N; ii++){
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk[jj*N+ii];   //weighted average
        }
    }
    
}

void combine_eig_coeffs( const fftw_complex *Jkx, const fftw_complex *Jky, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S){
    
    uint_t iwt; //index of weight to use
    fftw_complex Sk=0;
    for (uint_t ii=0; ii<N; ii++){
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            Sk=conj(Jkx[jj*N+ii])*Jky[jj*N+ii];
            
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk;   //weighted average
        }
    }
    
}

// NOTE: S is both an input AND output.
// Combines eigenspectra, fills S, but also computes the normalized absolute difference from what
// was previously in S.  This difference = || S_new - S_old ||/N is returned in diff.
void combine_eig_spec( const fftw_complex *Sk, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S, double &diff){
    
    uint_t iwt; //index of weight to use
    
    diff=0;
    double absdiff;
    fftw_complex Sl;

    for (uint_t ii=0; ii<N; ii++){
        Sl=S[ii];                               //stores previous value
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk[jj*N+ii];   //weighted average
        }
        absdiff=cabs(S[ii]-Sl);          //absolute difference
        diff=diff+absdiff;
    }
    
    diff=sqrt(diff)/N;
    
}

void combine_eig_coeffs( const fftw_complex *Jkx, const fftw_complex *Jky, uint_t K, uint_t N, const double *wk, uint_t nwk, fftw_complex *S, double &diff){
    
    uint_t iwt; //index of weight to use
    fftw_complex Sk;
    
    diff=0;
    double absdiff;
    fftw_complex Sl;
    for (uint_t ii=0; ii<N; ii++){
        Sl=S[ii];                               //stores previous value
        S[ii]=0;
        for ( uint_t jj=0; jj< K; jj++){
            Sk=conj(Jkx[jj*N+ii])*Jky[jj*N+ii];
            iwt=jj*nwk+(ii%nwk);
            S[ii]=S[ii]+wk[iwt]*Sk;   //weighted average
        }
        absdiff=cabs(S[ii]-Sl);       //absolute difference
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

