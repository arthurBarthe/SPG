function [ omega ] = Fourier_frequencies( N, delta )
%Returns the Fourier frequencies corresponding to a sample of N points with
%sampling step delta.
%Args:
%N                  int
%                   number of points
%delta              float
%                   sampling step
%Returns:
%omega              float[]
%                   Corresponding Fourier frequencies
omega=0:2*pi/delta/N:2*pi/delta*(1-1/N);
omega=fftshift(omega);
omega(1:floor(N/2))=omega(1:floor(N/2))-2*pi/delta;
end

