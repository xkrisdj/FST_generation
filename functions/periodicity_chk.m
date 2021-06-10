function [kx,ky,kz] = periodicity_chk(kx,ky,kz,Np,kk,dlx,dly,dlz,ifxp,ifyp,ifzp)
%PERIODICITY_CHK Summary of this function goes here
%   Detailed explanation goes here
%---------------------------------------------------------------------- 

% plot(ky,kz,'r.'); hold on
k2 = kk.^2;
seed=-randi([10 1000],20,20);
%     Add periodicity check
if (ifxp) 
    nmax = floor(kk*dlx/(2.0*pi));
    nmin = 1;
    if (nmax<nmin)
        disp(['Check allowed wavenumbers in FST'])
        disp(['nmax: ' num2str(nmax)])
        disp(['nmin: ' num2str(nmin)])        
        disp(['k   : ' num2str(kk)])
        error('Exiting')

    end  

for j=1:Np

    kn = sign(kx(j))*floor(abs(kx(j))*dlx/(2.0*pi));  % always make k smaller   

    if abs(kn)>nmax
        kn=kn-sign(kx(j));
    elseif abs(kn)==0
        kn=kn+sign(kx(j));
    end
    kx(j)=real(kn)*2.0*pi/dlx;

    flip = ran2(seed) ;           % coin toss
    if flip>0.5
        ky(j) = ky(j)  ;     % ky stays the same
        rtmp = k2-ky(j).^2.-kx(j).^2.;
    if rtmp>1
        kz(j) = sign(kz(j))*sqrt(rtmp);
    else
        rtmp = sqrt((k2-kx(j).^2.)/2.);
        ky(j) = sign(ky(j))*rtmp;
        kz(j) = sign(kz(j))*rtmp;
    end
    else
        kz(j) = kz(j);       % kz stays the same
        rtmp = k2-kx(j).^2.-kz(j).^2.;
    if rtmp>1                 
        ky(j) = sign(ky(j))*sqrt(rtmp);
    else
        rtmp = sqrt((k2-kx(j).^2.)/2.);
        ky(j) = sign(ky(j))*rtmp;
        kz(j) = sign(kz(j))*rtmp;
    end
    end       % flip
end        % j=1,Np

%----------------------------------

end
if (ifyp)
    nmax = fix(kk*dly/(2.0*pi)); %%% changed floor -> fix
    nmin = 1;
    if nmax<nmin
        disp(['Check allowed wavenumbers in FST'])
        disp(['nmax: ' num2str(nmax)])
        disp(['nmin: ' num2str(nmin)])        
        disp(['k   : ' num2str(kk)])
        error('Exiting')
    end

for j=1:Np                                          %%% changed floor -> fix
    kn = sign(ky(j))*fix(abs(ky(j))*dly/(2.0*pi));  % always make k smaller   

    if (abs(kn)>nmax)
        kn=kn-sign(ky(j));
    elseif abs(kn)==0
        kn=kn+sign(ky(j));
    end
    ky(j)=real(kn)*2.0*pi/dly;

    flip = ran2(seed(1,j));            % coin toss
    if flip>0.5
        kz(j) = kz(j);       % kz stays the same
        rtmp = k2-ky(j).^2.-kz(j).^2.;
        if rtmp>1
            kx(j) = sign(kx(j))*sqrt(rtmp);
        else
            rtmp = sqrt((k2-ky(j).^2.)/2.);
            kx(j) = sign(kx(j))*rtmp;
            kz(j) = sign(kz(j))*rtmp;
        end
    else
        kx(j) = kx(j);       % kx stays the same
        rtmp = k2-ky(j).^2.-kx(j).^2.;
        if rtmp>1      
            kz(j) = sign(kz(j))*sqrt(rtmp);
        else
            rtmp = sqrt((k2-ky(j).^2.)/2.);
            kx(j) = sign(kx(j))*rtmp;
            kz(j) = sign(kz(j))*rtmp;
        end
%     plot(ky(j),kz(j),'bo');
    end      % flip
end        % j=1,Np

%--------------------------

if (ifzp)
    nmax = fix(kk*dlz/(2.0*pi)); %%% changed floor -> fix
    nmin = 1;
    if nmax<nmin
        disp(['Check allowed wavenumbers in FST'])
        disp(['nmax: ' num2str(nmax)])
        disp(['nmin: ' num2str(nmin)])        
        disp(['k   : ' num2str(kk)])
        error('Exiting')
    end  

for j=1:Np                                          %%% changed floor -> fix
    kn = sign(kz(j))*fix(abs(kz(j))*dlz/(2.0*pi));  % always make k smaller   
    if (abs(kn)>nmax)
        kn=kn-sign(kz(j));
    elseif (abs(kn)==0)
        kn=kn+sign(kz(j));
    end
    kz(j)=real(kn)*2.0*pi/dlz;

    
    flip = ran2(seed(j,1));            % coin toss
    if (flip>0.5)
        kx(j) = kx(j);       % kx stays the same
        rtmp = k2-kx(j).^2.-kz(j).^2.;
        if rtmp>1
            ky(j) = sign(ky(j))*sqrt(rtmp);
        else
            rtmp = sqrt((k2-kz(j).^2.)/2.);
            kx(j) = sign(kx(j))*rtmp;
            ky(j) = sign(ky(j))*rtmp;
        end
    else
        ky(j) = ky(j);       % ky stays the same
        rtmp = k2-ky(j).^2.-kz(j).^2.;
    if rtmp>1      
        kx(j) = sign(kx(j))*sqrt(rtmp);
    else
        rtmp = sqrt((k2-kz(j).^2.)/2.);
        kx(j) = sign(kx(j))*rtmp;
        ky(j) = sign(ky(j))*rtmp;
    end
    end       % flip

end         % j=1,Np

end           % ifxp

end

