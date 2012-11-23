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

#include "dpss.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <sstream>

//usage:  dpss N nW [[lower_seq] upper_seq] [interp_method] [N_interp_base] [trace]
int main(int argc, char* argv[]);

//prints the workspace statistics using printf
void print_workspace(dpss_workspace work, char *prepend, char *append, char topchar, char bottomchar);
void print_workspace(dpss_workspace work);

template <class T> bool str_to(char *arg, T& dest);         //converts char array arguments to class T
void print_dpss(dpss *mydpss);                               //prints the results of dpss to stdout
void print_dpss(const double *h, uint_t n, const double *l, int nseq);  //prints the results of dpss to stdout
dpss_workspace cmd_parse_args(char *inargs[], int nargin);  //parses command line arguments, creating workspace

////////////////////////////////////////////////////////////////////////////
//  CODE
////////////////////////////////////////////////////////////////////////////

static const char helpstr[]={
"\n\
Usage: dpss n nW [[seql] sequ] [interp_method [interp_base]] [trace]\n\
\n\
`dpss' calculates and prints a set of discrete prolate spheroidal sequences\n\
(aka Slepian Sequences) and their corresponding concentrations of energy in\n\
the normalized bandwidth [-W, W].  These sequences are typically used in \n\
multitaper spectral analysis.\n\
    \n\
Examples:\n\
dpss 12 3 2        returns the first two DPSSs of length 12 and time-bandwidth\n\
                   product nW=3, along with their concentrations.\n\
dpss 12 3 3 5      returns the third through fifth DPSSs and concentrations.\n\
dpss 8388608 2.5 spline 32768   computes the first five DPSSs of length 2^15,\n\
                                nW=2.5, then interpolates to length 2^23 using\n\
                                natural cubic splines, re-normalizes and\n\
                                estimates the new energy concentrations.\n\
\n\
Required Parameters:\n\
  n     Length of the dpss sequences to be computed\n\
  nW    Time-bandwidth product for the sequences.  By definition, the first\n\
        dpss maximizes the concentration of energy in the normalized frequency\n\
        range f in [-W,W]. Typical values are nW = 2, 2.5, 3, 3.5 and 4.\n\
\n\
Optional Parameters:\n\
  sequ             The index of the upper dpss to be calculated.  Indices must\n\
                   fall in the range [1, N].  If `seql' is not defined, then\n\
                   this is the total number of sequences.  By default,\n\
                   sequ = floor(2nW). \n\
  seql             The index of the lower dpss to be calculated.  This must\n\
                   fall in the range [1, sequ].\n\
  interp_method    Interpolation method to employ when creating the sequences.\n\
                   By default, dpss will calculate sequences up to length\n\
                   n = 2^17-1 without interpolation.  For larger sequences,\n\
                   intermediate values must be interpolated.  Accepted methods\n\
                   are `linear' and `spline' (without quotes).  By default,\n\
                   interp_method=spline, which uses natural cubic splines.\n\
  interp_base      Length of the sequences from which to interpolate.  An\n\
                   interpolation method must be supplied to use this option.\n\
                   By default, interp_base = min( n, 2^17-1 ).\n\
  trace            trace = `trace' (without quotes) prints the method and \n\
                   computation parameters to the command window.\n\
  \n\
See the `mtpsd' Documentation for further details.\n\n"
};

//case insensitive comparison
int mystrnicmp (const std::string &s1, const std::string &s2, const int n)
{   
    for (int ii=0; ii<n; ii++){
        if (toupper(s1[ii])>toupper(s2[ii]) )
            return 1;
        else if ( toupper(s1[ii])<toupper(s2[ii]) )
            return -1;
    }
    return 0;
}

//converts string to type T using streams
template <class T>
bool str_to(char *arg, T& dest){
    std::istringstream sstr(arg);
    return !(sstr >> dest).fail();
}

//creates a string consisting of a repeated character (for display purposes)
char *rep_string(char *out, char rep, int n){
    
    if (n<0) n=0;       
    for (int ii=0; ii<n; ii++)
        out[ii]=rep;
    out[n]='\0';
    return out;
    
}

//prints workspace info
void print_workspace(dpss_workspace work, char *prepend, char *append, char topchar, char bottomchar){
    
    int  lmax=0,li=0,ln;
    char line[4][100];
    char tmp[100];
    
    sprintf(line[li++],"DPSS Summary ");
    sprintf(line[li++],"    Params: n=%d, nW=%.2f, Sequences dpss_%d to dpss_%d ", work.n, work.nW, work.seql+1, work.sequ+1);
    if (work.interp_method==NONE)
        sprintf(line[li++],"    Method: Tridiagonal formulation w/ even-odd splitting");
    else if (work.interp_method==LINEAR)
        sprintf(line[li++],"    Method: Interpolated from n=%d (Linear)",work.interp_base);
    else if (work.interp_method==SPLINE)
        sprintf(line[li++], "    Method: Interpolated from n=%d (Natural Cubic Splines)",work.interp_base);
    
    if (work.energy==true)
        sprintf(line[li++],"    Energy: Concentrations calculated using FFTs (FFTW)");

    for (int ii=0; ii<li; ii++)
        lmax=max<int>(lmax,strlen(line[ii]));
    
    int nprep, napp;
    nprep=strlen(prepend); 
    napp=strlen(append);
    
    //print top string
    if (topchar != '\0'){
        printf("\n%s\n",rep_string(tmp,topchar,lmax+nprep+napp));
        printf("%s%s%s\n",prepend,rep_string(tmp,' ',lmax),append);
    }
    
    //print all others
    for (int ii=0; ii<li; ii++){
        ln=strlen(line[ii]);
        printf("%s%s%s%s\n",prepend,line[ii],rep_string(tmp,' ',lmax-ln),append);    
    }
    
    //print bottom string
    if (bottomchar != '\0'){
        printf("%s%s%s\n",prepend,rep_string(tmp,' ',lmax),append);
        printf("%s\n\n",rep_string(tmp,bottomchar,lmax+nprep+napp));
    }
    
}

void print_workspace(dpss_workspace work){
    
    char *blank=(char*)"";
    print_workspace(work,blank,blank,'\0','\0');
    
}

//prints the discrete prolate spheroidal sequences and energies
void print_dpss(dpss *mydpss){
    
    print_dpss(mydpss->ph(), mydpss->size(0), mydpss->pl(), mydpss->size(1));
    
}

void print_dpss(const double *h, uint_t n, const double *l, int nseq){
    
    printf("dps_seq = [");
    for (unsigned int ii=0; ii<n; ii++){
        for (int jj=0; jj<nseq; jj++){
            printf(" %22.15e ",h[jj*n+ii]);
        }
        if (ii<n-1)
            printf(";\n\t   ");
        else
            printf("];\n\n");
    }
    
    printf("lambda = [ %22.15e", l[0]);
    for (int ii=1; ii<nseq; ii++){
        printf(" ;\n\t   %22.15e",l[ii]);
    }
    printf(" ];\n\n");
    
}

//parses the input arguments to create a workspace
//NOTE: little checking on valid parameters is done here
//      since the workspace will be fixed when a dpss class
//      is created.
dpss_workspace cmd_parse_args(char *inargs[], int nargin){
    
    dpss_workspace args;
    
    if (nargin < 2) throw ERR("Not enough arguments");
    
    uint_t n,tmp; double nW;
    
    int iarg=0;  int ii=0;    
    while ( (iarg <= INTERP_BASE) && (ii<nargin) ){
        
        switch(iarg){
            case SEQ_LENGTH:
                
                if (!str_to<uint_t>(inargs[ii],n))
                    throw ERR("Invalid dpss length");
                
                args.n=n;
                if ( n > MAX_N ){
                    args.interp_method=SPLINE;
                    args.interp_base=MAX_N;
                }           
                else {                 
                    args.interp_method=NONE;
                    args.interp_base=n;
                }
                ii++; iarg++;
                break;
                
            case TIME_HALF_BANDWIDTH:
                if (!str_to<double>(inargs[ii],nW))
                    throw ERR("Invalid nW");
                
                //defaults in case only 2 arguments supplied
                args.nW=nW;
                args.seql=0; args.sequ=floor(2*nW)-1;            
                ii++; iarg++;
                break;
                
            case SEQL:
                
                if (str_to<uint_t>(inargs[ii],tmp)){
                    args.seql=0;                //assume only upper bound given by default
                    args.sequ=max<uint_t>(tmp,1)-1;            //adjust so that numbers start at zero
                    ii++;
                }
                else 
                    iarg++;                    //pass to interpolation method
                
                iarg++;
                break;
                
            case SEQU:
                if (str_to<uint_t>(inargs[ii],tmp)){
                    args.seql=args.sequ;               //assume only one number given by default
                    args.sequ=max<uint_t>(tmp,1)-1;    //adjust so that numbers start at zero
                    ii++;
                }
                iarg++;
                break;
                
            case INTERP_METHOD:                 
                if ( mystrnicmp(inargs[ii],"linear", 6)==0 ){
                    args.interp_method=LINEAR;
                }
                else if ( (mystrnicmp(inargs[ii],"spline", 6)==0) || (mystrnicmp(inargs[ii],"interp", 6)==0))
                    args.interp_method=SPLINE;
                else {
                    throw ERR("Invalid interp method");
                }
                ii++; iarg++;
                break;
                
            case INTERP_BASE: 
                if (str_to<uint_t>(inargs[ii],tmp)){
                    if ( tmp < MAX_N)
                        args.interp_base = tmp;
                    else
                        args.interp_base=MAX_N;   
                    ii++;
                }
                iarg++;                
                break;
            default:
                throw ERR("Error parsing arguments");
        }
        
    }
    
    args.energy=true;   //always compute the energies in command line (not very expensive)
    
    return args;
    
}

//main routine
int main(int argc, char* argv[]){
    
    dpss_workspace dpss_args;
    bool trace=false;           //prints out computation stats
    
    if (argc>1){
        if ( ( mystrnicmp(argv[1],"--help", 6)==0 ) || ( mystrnicmp(argv[1],"?",1)==0 ) ){
            printf("%s",helpstr);
            return 0;
        }
    }

    if (mystrnicmp(argv[argc-1],"trace",5)==0){
        trace=true;
        argc=argc-1;
    }
    
    //parse arguments
    try{
        dpss_args=cmd_parse_args(&argv[1], argc-1);
    }
    catch(ERR err){
        printf("Error: %s\n", err.getmsg() );
		printf("%s",helpstr);
        return -1;
    }
    
    // Uses energize() to compute eigenvalues instead
    bool benergy=dpss_args.energy;
    dpss_args.energy=false;
    
    dpss myseq(dpss_args);
    try{
        myseq.compute();
    }
    catch(ERR err){
        printf("Error: %s\n", err.getmsg() );
        return -1;             
    }
    if (benergy==true)
        myseq.energize();
    
    
    if (trace==true){
        char *prep=(char*)"%  ";
        char *app=(char*)"  %";
        print_workspace(myseq.getinfo(),prep,app,'%','%');
    }
       
    print_dpss(&myseq);
    
    return 0;
    
}
