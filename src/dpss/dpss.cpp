/*

Copyright (c) 2010 C. Antonio Sanchez
    
This file is part of DPSS.
    
DPSS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
DPSS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
    
You should have received a copy of the GNU General Public License
along with DPSS.  If not, see <http://www.gnu.org/licenses/>.
    
*/  

// REFER TO dpss.h FOR A DESCRIPTION OF THE ROUTINES

#include "dpss.h"
#include <cmath>
#include <cstring>

//normalizes a vector h
void normalize_vec(double *h, uint_t n){
    
    double norm=0;
    for ( uint_t ii=0; ii<n; ii++)
        norm=norm+h[ii]*h[ii];
    norm=sqrt(norm);
    
    for ( uint_t ii=0; ii<n; ii++)
        h[ii]=h[ii]/norm;
    
}

//linearly interpolates to a new size 
void linear_interp( double *y, uint_t n, uint_t nout, double *z){
    
    double xi,xj,m,b;
    double *yin=new double[n];
    
    // make copy of input in case y=z;
    memcpy(yin,y,sizeof(double)*n);
    
    uint_t bin;
    for (uint_t jj=0; jj<nout; jj++){
        xj=(2.0*jj+1.0)/2.0/nout;
        
        if (xj<0.5/n) bin=0;
        else if (xj>1-0.5/n) bin=n-2;
        else bin=floor(n*xj-0.5);
        
        xi=(2.0*bin+1.0)/2.0/n;
        
        m=n*(y[bin+1]-y[bin]);
        b=y[bin]-m*xi;
        
        z[jj]=m*xj+b;
    }
    
    delete [] yin;
    
}

//interpolates to a new sizes using natural cubic splines
void spline_interp( double *y, uint_t n, uint_t nout, double *z ){
 
    double *c=new double[n];
    double *d=new double[n];
    
    //TRIDIAGONAL SOLVE
    //forward sweep
    c[0]=0.5; d[0]=0.5*3.0*n*(y[1]-y[0]);
    for ( uint_t ii=1; ii<(n-1); ii++ ){
        c[ii]= 1.0/(4.0-c[ii-1]);
        d[ii]= ( 3.0*n*(y[ii+1]-y[ii-1]) - d[ii-1] ) * c[ii];
    }
    c[n-1]=0;
    d[n-1]=( 3.0*n*( y[n-1]-y[n-2] ) - d[n-2])/(2.0-c[n-2]);
    
    //back substitution
    for (uint_t ii=n-2; ii>0; ii--){
        d[ii]=d[ii]-c[ii]*d[ii+1];
    }
    d[0]=d[0]-c[0]*d[1];    //last case separate because ii is uint_t (otherwise for loop goes on forever)
    
    // make copy of input in case y=z;
    memcpy(c,y,sizeof(double)*n);
    
    //INTERPOLATION USING PREVIOUSLY FOUND COEFFICIENTS
    //coefficients are stored in d
    //compute interpolation on uniform grid (actually extrapolates a little at the edges :S)
    double xi,xj,dy,dxn;
    uint_t bin;
    
    double aa,bb,cc,dd;                     //coefficients of cubic
    for ( uint_t jj=0; jj<nout; jj++ ){       
        xj=(2.0*jj+1.0)/2.0/nout;
        
        if (xj<0.5/n) bin=0;
        else if (xj>1-0.5/n) bin=n-2;
        else bin=floor(n*xj-0.5);
        xi=(2.0*bin+1.0)/2.0/n;
        
        dxn=(xj-xi)*n;
        dy=c[bin+1]-c[bin];
        
        aa=c[bin];
        bb=d[bin];
        cc= 3*dy*n-2*d[bin]-d[bin+1];
        dd=-2*dy*n+d[bin]+d[bin+1];       
        z[jj]=  aa + bb*dxn/n + cc*dxn*dxn/n + dd*dxn*dxn*dxn/n;        
    }
    
    delete [] c;
    delete [] d;
    
}

//polarizes the sequences
void polarize_dpss(double *h, uint_t n, uint_t iseq){
    
    double sum=0;
    
    //force positive mean for symmetric sequences, positive first lobe for antisymmetric
    int odd=(iseq%2);
    for (uint_t ii=0; ii<n/2; ii++)
        sum+=(odd*(n-1-2*ii)+1-odd)*h[ii];
        
    //if mean is negative, flip vector
    if (sum<0){
        for (uint_t ii=0; ii<n; ii++)
            h[ii]=-h[ii];
    }
}

//uses Relatively Robust Representations
void eig_rrr(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec){
    
    eig_rrr(n,D,E,il,iu,eig_val,eig_vec,n);
    
}

void eig_rrr(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec, uint_t vec_length){

    //output params
    char job='V', range='I', safe='S';   //find eigenvectors in index range using safe tolerance
    uint_t  m=iu-il+1;                   //number of returned eigenvalues/vectors
    uint_t ldz=vec_length;               //size of return array
    double vl,vu;                        //dummy variables    
    double abstol=2*dlamch_(&safe);      //safe tolerance
    const uint_t K=m;
    
    //outputs
    uint_t *ISUPPZ=new uint_t[2*K];     //support of Z
    uint_t info;                        //reports error status
    
    //workspace
    uint_t liwork=10*n, lwork=20*n;
    double *WORK = new double[lwork];
    uint_t *IWORK=new uint_t[liwork];
    
    //use LAPACK  function to find eigenvalues
    dstevr_(&job, &range, &n, D, E, &vl, &vu, &il, &iu, &abstol, &m, eig_val, eig_vec, &ldz, ISUPPZ, WORK, &lwork, IWORK, &liwork, &info );
    
    //cleanup
    delete [] WORK;  
    delete [] IWORK;
    delete [] ISUPPZ;
    
    if (info !=0 ){
        throw LAPACK_ERROR("DSTEVR failed");
    }
    
}

//uses inverse iteration (more expensive, but better accuracy)
void eig_iit(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec){
    
    eig_iit(n,D,E,il,iu,eig_val,eig_vec,n);
    
}

void eig_iit(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec, uint_t vec_length){
       
    //output params
    char range='I', order='B', safe='S';    //index range of eigenvalues, using block splitting and safe tolerance   
    double vl,vu;                           //dummy variables
    uint_t ldz=vec_length;                     //size of z array
    double abstol=2*dlamch_(&safe);         //safe minimum absolute tolerance
    const int K=iu-il+1;                 //number of returned eigenvalues/vectors
    
    //intermediaries
    uint_t m, nsplit;
    uint_t *IBLOCK=new uint_t[n];
    uint_t *ISPLIT=new uint_t[n];
    uint_t *IWORK=new uint_t[3*n];
    double *WORK = new double[5*n];               //dstebz needs 4n, but dstein needs 5n
    bool bdstebzfail=false;
    
    //outputs
    uint_t *IFAIL=new uint_t[K];
    uint_t info;
    
    //calculate eigenvalues
    dstebz_(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, D, E, &m, &nsplit, eig_val, IBLOCK, ISPLIT, WORK, IWORK, &info);
    if (info == 0)
        dstein_(&n, D, E, &m, eig_val, IBLOCK, ISPLIT, eig_vec, &ldz, WORK, IWORK, IFAIL, &info );
    else 
        bdstebzfail=true;
    
    //cleanup
    delete [] IBLOCK;  delete [] ISPLIT;
    delete [] IWORK;   delete [] WORK;
    delete [] IFAIL;
    
    if (info!=0){
        if (bdstebzfail==true)
            throw LAPACK_ERROR("DSTEBZ failed");
        throw LAPACK_ERROR("DSTEIN failed");
    }
}

void dpss_calc(uint_t n, double nW, int seql, int sequ, void(*eig_calc)(uint_t, double*, double*, uint_t, uint_t, double*, double*, uint_t), double *h){
    
    uint_t K=sequ-seql+1;
    bool is_n_even=(n%2==0);
    
    //sizes of even/odd problems
    uint_t n_even=ceil(n/2.0);
    uint_t n_odd=floor(n/2.0);  
    int seql_even=ceil(seql/2.0);
    int seql_odd=floor(seql/2.0);
    int sequ_even=floor(sequ/2.0);
    int sequ_odd=floor((sequ-1)/2.0);
    int nseq_even=sequ_even-seql_even+1;
    int nseq_odd=sequ_odd-seql_odd+1;
    
    //use same vectors for even/odd problem   
    double *D=new double[n_even];        
    double *E=new double[n_even];
    double *eig_val=new double[K];
    
    //fill tridiagonal matrix
    double cc=cos(2.0*PI*nW/n)/4.0;
    for (uint_t ii=0; ii < n_even; ii++){
        D[ii]=(n-1-2.0*ii)*(n-1-2.0*ii)*cc;
        E[ii]=(ii+1)*(n-ii-1)/2.0;
    }
    
    //central values for splitting
    double dc=D[n_even-1], ec=E[n_even-1];
    
    //EVEN PROBLEM
    if (is_n_even) 
        D[n_even-1]=dc+ec;
    else
        E[n_odd-1]=sqrt(2)*E[n_odd-1];
    
    uint_t start_odd=0,start_even=0;
    if (sequ%2==0)
        start_odd=n;
    else
        start_even=n;
    
    //calculate eigenvectors
    if (nseq_even>0)
        eig_calc(n_even,D,E,n_even-sequ_even,n_even-seql_even,eig_val,&h[start_even],2*n);
        
    //ODD PROBLEM
    D[n_even-1]=dc-ec;
    if (nseq_odd>0)
        eig_iit(n_odd,D,E,n_odd-sequ_odd,n_odd-seql_odd,&eig_val[nseq_even],&h[start_odd],2*n);
            
    //reverse eigenvector order
    double tmp;
    for (uint_t ii=0; ii<floor(K/2); ii++){
        for (uint_t jj=0; jj < n_even; jj++){
            tmp=h[(K-1-ii)*n+jj];
            h[(K-1-ii)*n+jj]=h[ii*n+jj];
            h[ii*n+jj]=tmp;
        }
    }
    
    //construct full vectors
    int c=1;
    if (seql%2==1) c=-1;
    for (uint_t ii=0; ii<K; ii++){
        for (uint_t jj=0; jj<n_odd; jj++){
            h[(ii+1)*n-1-jj]=c*h[ii*n+jj];
        }
        //correct centre value
        if (!is_n_even){
            if ((ii+seql)%2==1)
                h[ii*n+n_odd]=0;
            else
                h[ii*n+n_odd]=sqrt(2)*h[ii*n+n_odd];
        }
        c=-c;
    }
    
    //normalize and polarize the eigenvectors
    for (uint_t ii=0; ii<K; ii++){
        normalize_vec(&h[n*ii],n);
        polarize_dpss(&h[n*ii],n,seql+ii);
    }
    
    //cleanup
    delete [] D;  
    delete [] E;    
    delete [] eig_val;
}

//Reduces the problem using simple even/odd splitting (exploiting double symmetry)
void dpss_calc(uint_t n, double nW, int seql, int sequ, double *h){
     
    dpss_calc(n,nW,seql,sequ,&eig_iit, h);
    
}

// Adjusts settings in the dpss_workspace so that no information conflicts
// and the sequences and be calculated without throwing errors
void fix_workspace(dpss_workspace& work){
    
    if (work.n<0) work.n=0;
    
    if ( work.interp_base > MAX_N )
        work.interp_base = MAX_N;
    
    if ( ( work.n > MAX_N ) && ( work.interp_method == NONE ) ){
        work.interp_method=SPLINE;
        work.interp_base = min<uint_t>( MAX_N, work.interp_base );
    }

    //turn off (on) interpolation if not (is) required.
    if (work.interp_base == work.n)
        work.interp_method=NONE;
    else if ( work.interp_method == NONE )
        work.interp_method=SPLINE;   

    if (work.seql>work.sequ){
        uint_t tmp=work.seql;
        work.seql=work.sequ;
        work.sequ=tmp;
    }
    work.sequ=max<uint_t>(min<uint_t>(work.sequ,work.interp_base-1),0);
    work.seql=max<uint_t>(min<uint_t>(work.seql,work.interp_base-1),0);
    work.K=work.sequ-work.seql+1;
    
    //force nW to be between 1 and interp_base/2
    if (work.nW<1 || work.nW > work.interp_base/2.0) 
        work.nW=1;
    
}

//basic constructor
dpss::dpss(uint_t n, double nW){

    dpss_workspace *myinfo=&(this->info);

    myinfo->n=n;
    myinfo->nW=nW;
    myinfo->seql=0;
    myinfo->sequ=floor(2*nW)-2;         //2*nW-1 sequences by default
    myinfo->K=info.sequ-info.seql+1;
    myinfo->interp_method=NONE;
    myinfo->interp_base=n;
    myinfo->energy=false;
    
    fix_workspace(*myinfo);
   
    this->h = new double [(this->info.n)*(this->info.K)];
    this->l = new double [this->info.K];
        
}

//advanced constructor
dpss::dpss(dpss_workspace ws){

    this->info=ws;
    fix_workspace(this->info);
    this->h = new double [(this->info.n)*(this->info.K)];
    this->l = new double [this->info.K];
    
}

//Destructor
dpss::~dpss(){

    //cleanup filter and energy
    delete [] this->h;
    delete [] this->l;
    
}

//performs computation AND interpolation if required
void dpss::compute(){

    dpss_workspace *myinfo=&(this->info);
    double *hcmp;   //pointer to location sequences are computed
    bool barraysep=false;
    bool binterp=false;
    
    uint_t n=myinfo->n;
    uint_t K=myinfo->K;
    uint_t nb=myinfo->interp_base;
    uint_t istart=0;
    
    if (myinfo->interp_method != NONE) {
        
        if (nb > n) {
            barraysep=true;
            hcmp = new double[K*nb];
        }
        else{
            hcmp=this->h;
            istart=(n-nb)*K;
        }
        binterp=true;
    }
    else
        hcmp=this->h;

    //calculate sequences of length interp_base (if interp_method==NONE, then is same as n)    
    dpss_calc(nb, myinfo->nW, myinfo->seql,  myinfo->sequ,  &hcmp[istart]);
    
    //compute energies (or fake it by setting them all equal to unity)
    //This is done prior to interpolation because the computed sequences are
    //more reliable, particularly if linear interpolation is used.  To calculate
    //an estimate of the concentrations for the interpolated sequence (which
    //should be similar), run dpss::energize().
    if (myinfo->energy==true){
        compute_energy_concentrations(&hcmp[istart], nb, K, myinfo->nW, this->l);
    }
    else{
        for (uint_t ii=0; ii < myinfo->K; ii++)
            this->l[ii]=1;
    }
    
    //interpolate if required
    if (binterp==true){     
        void (*interp)(double*,uint_t,uint_t, double*);            
        if (myinfo->interp_method==LINEAR)
            interp=&linear_interp;
        else
            interp=&spline_interp;

        for (uint_t ii=0; ii<K; ii++){
            interp( &hcmp[istart+ii*nb], nb, n, &this->h[ii*n] );
            normalize_vec( &this->h[ii*n], n);
        }
    }

    //potential cleanup
    if (barraysep==true)
        delete [] hcmp;

}

//Computes energies of current sequences (possibly large due to interpolation)
void dpss::energize(){
    
    compute_energy_concentrations(this->h,this->info.n, this->info.K, this->info.nW,this->l);
    this->info.energy=true;
    
}

//length of sequence
uint_t dpss::length(){
    return (this->size(0))*(this->size(1));
}

//can return either length of sequence, or number of sequences
uint_t dpss::size(int dim){

    if (dim<1)
        return this->info.n;

    return this->info.K;

}

//index k lies in the range [seql, sequ]
double dpss::lambda(uint_t k){

    if ( ( k<= this->info.sequ) && (k >= this->info.seql) )
        return l[k-this->info.seql];
    return 0;
    
}

//returns the ith value of sequence k
//index k lies in the range [seql, sequ]
double dpss::operator()(uint_t k, uint_t i){

    if ( ( k<= this->info.sequ) && (k >= this->info.seql) && (i < this->info.n) && (i >= 0) )
        return h[ (k-this->info.seql)*(this->info.n) + i];
    return 0;
    
}
    
dpss_workspace dpss::getinfo(){
    return this->info;
}

//returns const pointer to filter h
const double* dpss::ph(){
    const double *mypoint = this->h;
    return mypoint;
}

//returns const pointer to energies l
const double* dpss::pl(){
    const double *mypoint = this->l;
    return mypoint;
}            
