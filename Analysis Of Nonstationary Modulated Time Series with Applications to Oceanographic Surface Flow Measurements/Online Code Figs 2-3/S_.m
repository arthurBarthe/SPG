function [ ESF ] = S_(ker, params, delta )
%This function returns the theoretical expected periodogram for the given
%parameters of the OU and Mat�rn, and the kernel passed to it.
%Parameters:
%ker [1xN]: kernel that corresponds to the changing coriolis frequencies.
%           Can be computed using the function coriolis_kernel.
%params [1x6]:  parameters for the OU+Mat�rn. Note that the second
%               parameter, params(2), is the shift for the coriolis 
%               frequencies (these are already accounted for in the kernel).
%               Ask me about this point if you don't understand.
%Output:
%ESF [1xN]: the expected periodogram at Fourier frequencies (from -pi to
%pi).
N = length(ker);
acv = maternacvs(params(4:6),N,delta) .* (1-(0:N-1)/N);
acv = acv + complexouacvs(params(1:3), N, delta) .* ker;
ESF2=2*fft(acv)-acv(1); ESF=abs(real(fftshift(ESF2)));
end

