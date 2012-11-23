%     Examples:
          Fs= 60;
          t = 0:1/Fs:5;
          x = cos(2*pi*20*t) + 0.2*randn(size(t));
          [S1,f,Sc1] = mtpsd( x,2.5,[],2*length(t),Fs);
          figure(1);
          plot( f, 10*log10([S1, Sc1]) );
%     Computes and plots the multitaper estimate using adaptive
%     weighting and the default 4 tapers.

          [S2,f,Sc2,F] = mtpsd( x,2.5,[],2*length(t),Fs,[],1-1/length(x));
          f_significant= f( F(:,1)>F(:,2) );
          printf('Significant frequencies:\t');
          disp(f_significant');
%     Computes and displays significant frequencies based on the F-test
%     at the 1% significance level.

          z = exp(I*2*pi*20*t) + 1/sqrt(2)*(1+I) + 0.2*( randn(size(t)) + I*randn(size(t)) );
          S3a = mtpsd( z,2.5,[],2*length(t),Fs, 'eigen', 'mean');
          S3b = mtpsd( z,2.5,[],2*length(t),Fs, 'eigen');
          figure(2);
          plot( f,10*log10([S3a S3b]));
%     Computes the spectrum of the complex series, z, weighting by
%     eigenvalues of tapers (energy concentrations). For S3a, the mean
%     is left in the data for spectrum calculations.  It is removed in
%     the computation of S3b.
