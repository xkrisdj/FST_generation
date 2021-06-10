function [kx, ky, kz]=period(kx1,ky1,kz1,Np,kk1,dlx,dly,dlz,ifxp,ifyp,ifzp)


if (ifyp && ifzp)

alpha_0 = 2*pi/dly;
ky0=2*pi/dly;
kz0=2*pi/dlz;

for j=1:Np
    
    ky(j) = floor(ky1(j)/ky0+0.5)*ky0;
    kz(j) = floor(kz1(j)/kz0+0.5)*kz0;
    
%     rtmp = kk1-ky(j).^2-kz(j).^2;     
%     kx(j) =sign(kx1(j))*(sqrt(rtmp));
    kx(j) =kx1(j);
end


else
    for j=1:Np
    kx(j)=kx1(j);
    ky(j)=ky1(j);
    kz(j)=kz1(j);
    end
end %end if periodic