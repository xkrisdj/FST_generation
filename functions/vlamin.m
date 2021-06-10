function [VLAMIN] = vlamin(VEC,N)
%VLAMIN Summary of this function goes here
%   Detailed explanation goes here

TMIN = 99.0E20;

for I=1:N
TMIN = min(TMIN,abs(VEC(I)));
end

VLAMIN = TMIN;
end
