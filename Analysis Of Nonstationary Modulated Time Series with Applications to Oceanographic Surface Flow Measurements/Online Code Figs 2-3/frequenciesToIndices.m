function [ index ] = frequenciesToIndices(freq, omega)
%Converts a frequency with a given unit to index.
%Args:
%freq                   float
%                       A frequency to convert to indices.
%omega                  Fourier frequencies. Unit should be the same as
%                       that of freqs.
%We treat positive and negative frequencies differently.
temp = find(omega > freq);
if freq > 0
    index = temp(1)-1;
else
    index = temp(1)+1;
end
end

