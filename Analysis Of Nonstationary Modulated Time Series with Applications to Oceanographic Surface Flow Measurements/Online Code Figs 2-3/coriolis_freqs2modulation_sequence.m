function [ g ] = coriolis_freqs2modulation_sequence( freqs, delta )
%Returns the modulation sequence that corresponds to a sample of coriolis
%frequencies, given the sampling step delta. The units of the frequencies
%and of delta must correspond. For instance if freqs is in cycles per day,
%delta should be in days.
%Args:
%freqs                  float[]
%                       Time series of coriolis frequencies
%delta                  float[]
%                       The sampling step of the time series
g = exp(1i*cumsum(freqs .* delta));
end

