function [ freqs ] = coriolis_frequency(latitude, unit)
%This function returns the Coriolis frequency that corresponds to the
%passed latitude(s). Optional argument unit specifies what unit the
%frequency should be in.
%Args:
%latitude                   float[]
%                           latitude or array of latitudes in degrees
%unit (optional)            string
%                           A string specifying the desired unit. Permitted
%                           values are 'r-h' (radians per hour, default), 
%                           'r-d'(radians per day), 'c-h' (cycles per hour),
%                           'c-d' (cycles per day)
%Returns:
%freqs                      float[]
%                           Corresponding Coriolis frequencies
switch nargin
    case 1
        unit = 'r-h';
end
latitude_radians = latitude / 180 * pi;
freq_r_h = -2*7.2921*10^(-5)*3600*sin(latitude_radians);
switch unit
    case 'r-h'
        freqs = freq_r_h;
    case 'r-d'
        freqs = 24 * freq_r_h;
    case 'c-d'
        freqs = freq_r_h * 24 / 2 / pi;
    case 'c-h'
        freqs = freq_r_h / 2 / pi;
end
end

