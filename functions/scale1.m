function [xs,ys,zs] = scale1(x,y,z,Np,r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% xs=zeros(Np,1);
% ys=zeros(Np,1);
% zs=zeros(Np,1);

for i=1:Np
xs(i)=x(i)*r;
ys(i)=y(i)*r;
zs(i)=z(i)*r;
end


end

