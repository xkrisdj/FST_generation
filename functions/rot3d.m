function [xr,yr,zr] = rot3d(x,y,z,Np,rotx,roty,rotz)
%ROT3D Summary of this function goes here
%   Detailed explanation goes here

xr=zeros(Np,1);
yr=zeros(Np,1);
zr=zeros(Np,1);

cx=cos(rotx);
sx=sin(rotx);
cy=cos(roty);
sy=sin(roty);
cz=cos(rotz);
sz=sin(rotz);

for i=1:Np
xx=x(i);
yy=y(i);
zz=z(i);

xr(i)=xx*cy*cz+yy*cy*sz-zz*sy;
yr(i)=xx*(cz*sx*sy-cx*sz)+yy*(cx*cz+sx*sy*sz)+zz*cy*sx;
zr(i)=xx*(cx*cz*sy+sx*sz)-yy*(cz*sx-cx*sy*sz)+zz*cx*cy;

end

end

