function [ freq_ ] = convertFrequency(freq, unit, unit_)
%Returns the value of the frequency converted into a different unit.
%Args:
%freq                   float
%                       frequency value
%unit                   str
%                       Unit of the frequency value provided by freq.
%                       Authorized values are of the form 'A per B' where A
%                       can take values 'radians' or 'cycles', B can take
%                       values 'second', 'minute', 'hour', 'day'.
%unit_                  str
%                       Unit to convert to. Authorized values are the same
%                       as for unit.
%Returns
%freq_                  float
%                       The frequency converted in the unit specified by
%                       argument unit_.
t = find(unit == ' ');
unit1 = unit(1:t(1)-1);
unit2 = unit(t(2)+1:end);
t = find(unit_ == ' ');
unit1_ = unit_(1:t(1)-1);
unit2_ = unit_(t(2)+1:end);
k = {'radians', 'cycles'};
v = {2*pi, 1};
k2 = {'second', 'minute', 'hour', 'day'};
v2 = {1, 60, 3600, 3600*24};
conversion = containers.Map(k,v);
conversion2 = containers.Map(k2, v2);
freq_ = freq * conversion2(unit2_) / conversion2(unit2);
freq_ = freq_ * conversion(unit1_) / conversion(unit1);
end

