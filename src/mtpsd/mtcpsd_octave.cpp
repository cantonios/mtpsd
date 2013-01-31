/*

Copyright (c) 2010 C. Antonio Sanchez

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

#include <octave/oct.h>
#include <octave/CMatrix.h>

#include "mtcpsd.h"
#include <cmath>

//Order of arguments:  x, nW,     ,nseq,   nfft, fs,  p_conf, p_ftest, 'mean', 'method'
//                 or     tapers, lambda                           
//
// Possible outputs (in order);  S, f, Sc, F, Jk, wk, tapers, lambda

const char helpstr[]= {\
"-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW})  \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW}, @var{nseq})  \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW}, @var{nseq}, @var{NFFT}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW}, @var{nseq}, @var{NFFT}, @var{Fs}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW}, @var{nseq}, @var{NFFT}, @var{Fs}, @var{pc}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{nW}, @var{nseq}, @var{NFFT}, @var{Fs}, @var{pc}, @var{pf}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{h}, @var{l}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{h}, @var{l}, @var{NFFT}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{h}, @var{l}, @var{NFFT}, @var{Fs}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{h}, @var{l}, @var{NFFT}, @var{Fs}, @var{pc}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@var{x}, @var{y}, @var{h}, @var{l}, @var{NFFT}, @var{Fs}, @var{pc}, @var{pf}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@dots{}, @var{'method'}) \n\
@deftypefnx {Loadable Function} {@var{S} =} mtcpsd (@dots{}, @var{'mean'}) \n\
@deftypefnx {Loadable Function} {[@var{S}, @var{f}, @var{Sc}, @var{Jkx}, @var{Jky}, @var{wk}, @var{h}, @var{l}] =} mtcpsd (@dots{}) \n\
@cindex dpss slepian multitaper psd \n\
\n\
Uses Thomson's Multitaper (MT) method to estimate the Cross-Power Spectral Density (CPSD) of two one-dimensional time-series. The method uses orthogonal tapers to obtain a set of uncorrelated spectrum estimates, called eigenspectra.  These are then combined either linearly (equal or eigenvalue weighting), or non-linearly (adaptive weighting) to form a single spectrum. The resulting estimate has been shown to have lower variance and broad-band bias properties than the classical periodogram [1]. \n\
\n\
The data tapers are taken to be discrete prolate spheroidal sequences because of their desirable frequency characteristics. These sequences can be specified by their time-bandwidth product, @var{nW}, or can be supplied as column vectors of @var{h} with corresponding energy concentrations @var{l}.  The larger the value of @var{nW}, the larger the main lobe of the spectral windows.  See the @code{dpss} documentation for more details regarding their definition and calculation.\n\
\n\
A particular weighting method for the multitaper estimate can be specified by the @var{'method'} string, which can take the following values:\n\
@table @asis\n\
@item @var{'equal'}\n\
Each eigenspectrum is weighted equally.\n\
@item @var{'eigen'}\n\
Each eigenspectrum is weighted by their energy concentration (eigenvalue).\n\
@end table\n\
By default, @code{mtcpsd} will use eigenvalue weighting.\n\
\n\
A confidence interval is computed if @var{Sc} is requested.  This interval is based on the assumption that the spectrum estimate follows a scaled chi-squared distribution, where the degrees of freedom is dependent on the eigenspectrum weights.  The width of the interval can be specified by the probability value @code{pc}.\n\
\n\
In order to reduce bias effects from constant terms, weighted means of the data are removed prior to computing each eigenspectrum.  These are never re-introduced, leaving estimation of the mean up to the user.  This behaviour can be overridden by including @var{'mean'} as the final input variable.\n\
\n\
\n\
Input Variables:\n\
\n\
@table @asis\n\
@item @var{x}, @var{y}\n\
One-dimensional time-series data (real or complex) of which to compute the cross power spectral density.\n\
@item @var{nW}\n\
The time-bandwidth product for the DPSSs.  Typical values are 2, 2.5, 3, 3.5 and 4.\n\
@item @var{nseq}\n\
The number of tapers to use.  The default value is \n\
@iftex\n\
@tex\n\
$\\lfloor 2nW \\rfloor - 1$,\n\
@end tex\n\
@end iftex\n\
@ifinfo\n\
floor(2@code{nW}),\n\
@end ifinfo\n\
since this is the number of DPSSs with energy concentrations close to one.\n\
@item @var{NFFT}\n\
The length of the FFT used in spectrum calculations.  The returned spectrum will consist of @var{NFFT} values at equally spaced frequencies.  The default value is @code{length}(@var{x}).\n\
@item @var{Fs}\n\
The sampling frequency.  The returned spectrum is evaluated in the Nyquist range [-Fs/2, Fs/2], beginning with the zero-frequency.  The default value is @var{Fs}=1.\n\
@item @var{pc}\n\
The width of the confidence interval.  The output @var{Sc} corresponds to the lower and upper limits of the \n\
@iftex\n\
@tex\n\
$pc\\times100$\%\n\
@end tex\n\
@end iftex\n\
@ifinfo\n\
@var{pc}*100%\n\
@end ifinfo\n\
confidence interval, assuming @var{S} follows a scaled chi-squared distribution.  The default value is @var{pc}=0.95.\n\
@item @var{h}\n\
User-supplied data tapers.  These are typically the first few discrete prolate spheroidal sequences, and can be computed using the @code{dpss} function.  They must have the same length as the time-series, @var{x}.  Tapers are taken to be the columns of @var{h}.\n\
@item @var{l}\n\
Vector of energy concentrations for the user-supplied tapers.  This must have @var{nseq} entries, where @var{nseq} is the number of columns in @var{h}.  By definition, @var{l}(i) is the ratio of power in [-W, W] to total power for the i-th taper, where W is a normalized frequency defined by @var{nW}/@code{length}(@var{x}). \n\
@item @var{'method'}\n\
The weight method used to combine the individual eigenspectra into a single estimate.  Valid entries are: @var{'equal'}, @var{'eigen'}, and @var{'adapt'}.\n\
@item @var{'mean'}\n\
If specified, leaves the mean value in the data when computing the spectrum.  Otherwise, weighted means are removed (weighted by the tapers) prior to computing each eigenspectrum.  This is done to force the zero-frequency component to zero, eliminating any bias caused by constant terms.\n\
@end table\n\
\n\
If the empty vector, @code{[]}, is used for @var{nseq}, @var{NFFT}, @var{Fs}, @var{pc}, or @var{pf}, the default value is used.\n\
\n\
\n\
Output Variables:\n\
\n\
@table @asis\n\
@item @var{S}\n\
The two-sided power spectral density estimate.  @var{S}(1) corresponds to the zero-frequency value.\n\
@item @var{f}\n\
Vector of frequencies.  @var{f} is in the range [0, @var{Fs}].\n\
@item @var{Sc}\n\
Confidence interval.  The first and second columns are the lower and upper bounds of the \n\
@iftex\n\
@tex\n\
$pc\\times100$\%\n\
@end tex\n\
@end iftex\n\
@ifinfo\n\
@var{pc}*100%\n\
@end ifinfo\n\
confidence interval, respectively.  It is assumed that @var{S} follows a scaled chi-squared distribution.\n\
@item @var{Jkx}, @var{Jky}\n\
The eigencoefficients.  These are formed by tapering the data and applying the FFT.\n\
@item @var{wk}\n\
The weight vectors.  @var{wk} has dimensions @var{NFFT}\n\
@iftex\n\
@tex\n\
$\\times$\n\
@end tex\n\
@end iftex\n\
@ifinfo\n\
x\n\
@end ifinfo\n\
@var{nseq}.  The i-th row is the set of weights used to calculate the i-th value of the spectrum: \n\
@iftex\n\
@tex\n\
@var{S}$(i)=\\sum_k$@var{wk}$(i)Sk$.\n\
@end tex\n\
@end iftex\n\
@ifinfo\n\
@var{S}=sum(@var{wk}.*Sk, 2).\n\
@end ifinfo\n\
@item @var{h}\n\
The data tapers, as columns.\n\
@item @var{l}\n\
Row vector containing the energy concentrations of the tapers.\n\
@end table\n\
\n\
REFERENCES\n\
\n\
[1] Percival, D.B., and A.T. Walden, Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques, Cambridge University Press, 1993.\n\
\n\
[2] Sanchez, C.A., @code{mtcpsd} Documentation, http://sourceforge.net/projects/mtcpsd, 2013.\n\
@end deftypefn\n"
};

enum PARSE_ARGS { DATA_LENGTH, TAPER_INFO, NUM_TAPERS, FFT_LENGTH, SAMPLE_FREQ, CONF_BOUND };

struct mtcpsd_args{
    mtcpsd_workspace info;
    bool compute_tapers;
    double p_conf;
};

mtcpsd_args parse_args(octave_value_list args, int nargout);

template <class T>
octave_value_list compute_spectrum(T *x, T *y, mtcpsd_args inargs, const double *h, const double *l, int nargout);
template  octave_value_list compute_spectrum<double>(double *x, double *y, mtcpsd_args inargs, const double *h, const double *l, int nargout);
template  octave_value_list compute_spectrum<fftw_complex>(fftw_complex *x, fftw_complex *y, mtcpsd_args inargs, const double *h, const double *l, int nargout);

//case insensitive comparison
int mystrnicmp(const std::string &s1, const std::string &s2, const int n)
{   
    for (int ii=0; ii<n; ii++){
        if (toupper(s1[ii])>toupper(s2[ii]) )
            return 1;
        else if ( toupper(s1[ii])<toupper(s2[ii]) )
            return -1;
    }
    return 0;
}

mtcpsd_args parse_args(octave_value_list args, int nargout){
    
    //Order of arguments:  x, nW,     [nseq]   optional:  nfft, fs,  p_conf, 'mean', 'method'
    //                 or     tapers, lambda                           
    //                                  
    
    int nargin = args.length();
    mtcpsd_workspace info;
    mtcpsd_args retval;
    
    if (nargin<2)
        throw ERR("Not enough input arguments");
    
    // Check data type and size
    if (args(0).is_complex_matrix()){
        info.dtype=COMPLEX_DATA;
        info.n=args(0).matrix_value().nelem();
    }
    else if (args(0).is_real_matrix()){
        info.dtype=REAL_DATA;
        info.n=args(0).matrix_value().nelem();
    }
    else
        throw ERR("Invalid data type");
    
    // fill in defaults
    info.N=info.n;
    info.Fs=1;
    info.nwk=1;                      // length of degree of freedom vector (1 or n)
    info.remove_mean=true;
    retval.p_conf=0.95;
    retval.compute_tapers=true;
    
    //check to see if last argument is 'mean':
    if (args(nargin-1).is_string()){
        if (mystrnicmp( args(nargin-1).string_value(), "mean", 4)==0 ) {
            info.remove_mean=false;
            nargin--;
        }
    }
    
    //checks if last (or second-last if 'mean' was given) argument supplies a method
    if (args(nargin-1).is_string()){
        if ( mystrnicmp( args(nargin-1).string_value(), "equal", 5)==0 )
            info.weight_method=EQUAL;
        else if ( mystrnicmp( args(nargin-1).string_value(), "eigen", 5)==0 )
            info.weight_method=EIGEN;
        else
            throw ERR("Invalid weighting method");
        
        nargin--;   //effectively remove entry
    }
    else
        info.weight_method=EIGEN;
    
    if (nargin<3)
        throw ERR("Not enough input arguments");
    
    int iarg=TAPER_INFO;
    int ii=2;
    
    uint_t utmp;    //temp variables
    double tmp;
    
    while ( (iarg <= CONF_BOUND) && (ii<nargin) ){
        
        switch(iarg){
            case TAPER_INFO:
                
                //check if tapers are supplied
                if (args(ii).is_matrix_type()){
                    
                    if ((uint_t)args(ii).matrix_value().dims()(0)!=info.n)
                        throw ERR("Invalid taper length");
                    
                    info.K=args(ii).matrix_value().dims()(1);                    
                    info.nW=-1;      //filler value;
                    ii++;
                    iarg++;
                    
                    if (ii<nargin){
                        if (args(ii).is_matrix_type()){
                            if ( (uint_t)args(ii).matrix_value().nelem() != info.K)
                                throw ERR("Energy vector has wrong length");
                            else {
                                ii++;
                                retval.compute_tapers=false;
                            }
                        }
                        else
                            throw ERR("Vector of energies required");
                            
                    }
                    else
                        throw ERR("Vector of energies required");
                    
                    
                }
                else if (args(ii).is_scalar_type()){
                    info.nW=args(ii).double_value();
                    info.K=floor(info.nW*2)-1;
                    ii++;
                }
                else
                    throw ERR("Invalid taper definition");
                
                iarg++;
                
                break;
            case NUM_TAPERS:
                if (args(ii).is_empty()){
                    //default
                    info.K=floor(info.nW*2)-1;
                    iarg++;
                    ii++;
                    break;
                }
                    
                if (!args(ii).is_scalar_type())
                    throw ERR("Invalid input.  Expected: N_TAPERS or NFFT");
                
                utmp=args(ii).int_value();
                
                if (utmp<info.N){
                    info.K=utmp;
                    ii++;
                }
                iarg++;
                
                break;
            case FFT_LENGTH:
                if (args(ii).is_empty()){
                    //default
                    info.N=info.n;
                    iarg++;
                    ii++;
                    break;
                }
                
                if (!args(ii).is_scalar_type())
                    throw ERR("Invalid FFT length");
                
                utmp=args(ii).int_value();
                if (utmp>=info.N){
                    info.N=utmp;
                    ii++;
                }
                iarg++;
                
                break;
            case SAMPLE_FREQ:
                if (args(ii).is_empty()){
                    //default
                    info.Fs=1;
                    iarg++;
                    ii++;
                    break;
                }
                
                if (!args(ii).is_scalar_type())
                    throw ERR("Invalid sample frequency");
                
                info.Fs=args(ii).double_value();
                ii++;
                iarg++;
                
                break;
            case CONF_BOUND:
                if (args(ii).is_empty()){
                    //default
                    retval.p_conf=0.95;
                    iarg++;
                    ii++;
                    break;
                }
                
                if (!args(ii).is_scalar_type())
                    throw ERR("Invalid confidence interval");
                
                tmp=args(ii).double_value();
                
                if ( (tmp <= 0 ) || ( tmp >= 100 ) )
                    retval.p_conf=0;
                else if ( tmp >= 1 )
                    retval.p_conf=tmp/100;
                else
                    retval.p_conf=tmp;
                ii++;
                iarg++;
                
                break;
            default:
                throw ERR("Invalid input arguments");
        }
    }
    fix_workspace(info);
    retval.info=info;
    
    return retval;
    
}

template <class T>
octave_value_list compute_spectrum(T *x, T *y, mtcpsd_args inargs, const double *h, const double *l, int nargout){

    octave_value_list retval;
    mtcpsd_workspace mt_work=inargs.info;
    
    mtcpsd<T> spec(x, y, h, l, mt_work);
    try{
        spec.compute();
    }
    catch(ERR err){
        printf("mtcpsd Error: %s\n", err.getmsg() );
        return retval;
    }
    
    // Build and return outputs
    // Possible outputs (in order);  S, f, Sc, Jkx, Jky, wk, tapers, lambda
    
    int nextarg=7;
    if (nargout > nextarg){
        Matrix lambdaout(1, mt_work.K);
        for ( uint_t ii=0; ii< mt_work.K; ii++){
            lambdaout(0,ii)=spec.lambda(ii);
        }
        retval(nextarg)=lambdaout;
    }
    nextarg--;
    if (nargout > nextarg){
        Matrix taperout(mt_work.n, mt_work.K);
        for ( uint_t ii=0; ii< mt_work.K; ii++){
            for (uint_t jj=0; jj<mt_work.n; jj++){
                taperout(jj,ii)=spec.taper(ii,jj);
            }
        }
        retval(nextarg)=taperout;
    }
    nextarg--;
    if (nargout > nextarg){
        Matrix wkout(mt_work.N,mt_work.K);
        for ( uint_t ii=0; ii< mt_work.K; ii++){
            
            double wt=spec.wt(ii,0);
            for (uint_t jj=0; jj<mt_work.N; jj++){
                wkout(jj,ii)=wt;
            }
        }
        retval(nextarg)=wkout;
    }
    nextarg--;
    if (nargout > nextarg){
        ComplexMatrix Jkxout(mt_work.N,mt_work.K);
        for ( uint_t ii=0; ii< mt_work.K; ii++){
            for (uint_t jj=0; jj<mt_work.N; jj++){
                Jkxout(jj,ii)=spec.eig_coeff_x(ii,jj);
            }
        }
        retval(nextarg)=Jkxout;
    }
    nextarg--;
    if (nargout > nextarg){
		ComplexMatrix Jkyout(mt_work.N,mt_work.K);
		for ( uint_t ii=0; ii< mt_work.K; ii++){
			for (uint_t jj=0; jj<mt_work.N; jj++){
				Jkyout(jj,ii)=spec.eig_coeff_y(ii,jj);
			}
		}
		retval(nextarg)=Jkyout;
	}
    nextarg--;
    if (nargout > nextarg){
        Matrix Scout(mt_work.N,2);
        double lfactor=spec.conf_factor(0,LOWER,inargs.p_conf);
        double ufactor=spec.conf_factor(0,UPPER,inargs.p_conf);
        
        for (uint_t jj=0; jj<mt_work.N; jj++){
            Scout(jj,0)=cabs(spec(jj))*lfactor;
            Scout(jj,1)=cabs(spec(jj))*ufactor;
        }
        retval(nextarg)=Scout;
    }
    nextarg--;
    if (nargout > nextarg){
        Matrix fout(mt_work.N,1);
        for (uint_t ii=0; ii<mt_work.N; ii++)
            fout(ii)=spec.freq(ii);
        retval(nextarg)=fout;
        
    }
    nextarg--;
    ComplexMatrix Sout(mt_work.N,1);
    for (uint_t ii=0; ii<mt_work.N; ii++)
        Sout(ii)=spec(ii);
    retval(0)=Sout;
    
    return retval;
    
}

//main interface routine
DEFUN_DLD (mtcpsd, args, nargout, helpstr)
{
    
    octave_value_list retval;
    
    mtcpsd_args inargs;
    mtcpsd_workspace mt_work;
    try{
        inargs = parse_args(args, nargout);
        mt_work=inargs.info;
    }
    catch (ERR err){
        printf("Error: %s\n", err.getmsg() );
        return retval;
    }
    catch(...){
        printf("Unknown error\n");
        return retval;
    }
    
    if (inargs.compute_tapers==true) {
        dpss_workspace taper_info;
        taper_info.n=mt_work.n;
        taper_info.nW=mt_work.nW;
        taper_info.seql=0;
        taper_info.sequ=mt_work.K-1;
        taper_info.K=mt_work.K;
        taper_info.interp_method=NONE;
        taper_info.interp_base=mt_work.n;
        taper_info.energy=true;
        fix_workspace(taper_info);
            
        dpss tapers(taper_info);
        try{
            tapers.compute();
        }
        catch(ERR err){
            printf("DPSS Error: %s\n", err.getmsg() );
            return retval;
        }
        const double *h=tapers.ph();
        const double *l=tapers.pl();
        if ( mt_work.dtype==REAL_DATA){
            retval=compute_spectrum<double>((double *)args(0).matrix_value().fortran_vec(), (double *)args(1).matrix_value().fortran_vec(),inargs, h, l, nargout);
        }
        else {
            retval=compute_spectrum<fftw_complex>( (fftw_complex *)args(0).complex_matrix_value().fortran_vec(), (fftw_complex *)args(1).complex_matrix_value().fortran_vec(), inargs, h, l, nargout);
        }        

    }
    else {
        
        //for some reason, repeating the const double *h=... here causes an error
        // I suspect that the *fortran_vec is deleted for some reason before the call to compute_spectrum
        // if declared out here
        if ( mt_work.dtype==REAL_DATA){
            retval=compute_spectrum<double>((double *)args(0).matrix_value().fortran_vec(),(double *)args(1).matrix_value().fortran_vec(), inargs, (double *)args(2).matrix_value().fortran_vec(), (double *)args(3).matrix_value().fortran_vec(), nargout);
        }
        else {
            retval=compute_spectrum<fftw_complex>( (fftw_complex *)args(0).complex_matrix_value().fortran_vec(), (fftw_complex *)args(1).matrix_value().fortran_vec(), inargs, (double *)args(2).matrix_value().fortran_vec(), (double *)args(3).matrix_value().fortran_vec(), nargout);
        }
    }
    
    return retval;
    
}
