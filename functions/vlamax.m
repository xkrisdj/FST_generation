function [VLAMAX] = vlamax(VEC,N)
%VLAMAX Summary of this function goes here
%   Detailed explanation goes here

TMAX = -99.0E20;

for I=1:N
TMAX = max(TMAX,abs(VEC(I)));
end

VLAMAX = TMAX;
end
