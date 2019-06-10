function [f,mx1,phase] = PowerSpectrum(XType1)

XType1 = XType1(:);

TimeFactor = numel(XType1);

XType1Normalized = XType1;

XType1NormalizedMean = mean(XType1Normalized,2);

sampling_freq = length(XType1(:,1));
Fs = sampling_freq; 

% Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
nfft= 2^(nextpow2(length(XType1NormalizedMean))); 

% Take fft, padding with zeros so that length(fftx) is equal to nfft 

fftx1 = fft(XType1NormalizedMean,nfft); 

% Calculate the numberof unique points
NumUniquePts = ceil((nfft+1)/2); 

% FFT is symmetric, throw away second half 

fftx1 = fftx1(1:NumUniquePts); 

phase = angle(fftx1);

mx1 = fftx1.*conj(fftx1)/length(fftx1)^2;

% This is an evenly spaced frequency vector with NumUniquePts points. 
numOfDays = length(XType1(:,1));
f = (0:NumUniquePts-1)*(Fs)/nfft; 

end