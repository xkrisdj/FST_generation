function [xt,yt,zt] = trans(x,y,z,Np,xx,yy,zz)
%TRANS Summary of this function goes here
%   Detailed explanation goes here
xt=zeros(Np,1);
yt=zeros(Np,1);
zt=zeros(Np,1);

for i=1:Np
xt(i)=x(i)+xx;
yt(i)=y(i)+yy;
zt(i)=z(i)+zz;
end

end

