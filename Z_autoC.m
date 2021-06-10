% clear all
% close all
% clc
% 
% load fieldC24u.mat
% ux=struct('U',{U.u});
% ux=cell2mat(struct2cell(ux));
% clear U
% 
% load fieldC24v.mat
% uy=struct('Uv',{Uv.v});
% uy=cell2mat(struct2cell(uy));
% clear Uv
% 
% load fieldC24w.mat
% uzz=struct('Uw',{Uw.w});
% uzz=cell2mat(struct2cell(uzz));
% clear Uw


% y=squeeze(y(1,1:end,1:end));
% z=squeeze(z(1,1:end,1:end));

y_y=linspace(0,dly,ny);
z_z=linspace(0,dlz,nz);
%% 

uyz=ux;

%% uy
Ny=size(uyz,1);
Nlag=Ny-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3)-0000 % loop over different times  20, 40, 60, 80, 90, 96, 98
    mm=mm+1;
    cc=0;
for m=1:size(uyz,2)
    cc = cc + 1 ;
    uu=uyz(:,m,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(:,m,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


uym=uzm;
uym_fft=uzmt_FFT;

disp(['Integral length scale uy: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,uym,'c-'); hold on
% plot(y_y,uym_fft,'c--')


%

%% uz
Nz=size(uyz,2);
Nlag=Nz-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3)-0000 % loop over different times 20, 40, 60, 80, 90, 96, 98
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=uyz(m,:,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(m,:,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


uzmm=uzm;
uzm_fft=uzmt_FFT;

disp(['Integral length scale uz: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,uzmm,'m-'); hold on
% plot(y_y,uzm_fft,'m--')


%% 

uyz=uy;

%% vy
Ny=size(uyz,1);
Nlag=Ny-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,2)
    cc = cc + 1 ;
    uu=uyz(:,m,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(:,m,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


vym=uzm;
vym_fft=uzmt_FFT;

disp(['Integral length scale vy: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,vym,'k-'); hold on
% plot(y_y,vym_fft,'k--')


%

%% vz
Nz=size(uyz,2);
Nlag=Nz-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=uyz(m,:,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(m,:,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


vzmm=uzm;
vzm_fft=uzmt_FFT;

disp(['Integral length scale vz: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,vzmm,'g-'); hold on
% plot(y_y,vzm_fft,'g--')


%% 

uyz=uzz;

%% wy
Ny=size(uyz,1);
Nlag=Ny-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,2)
    cc = cc + 1 ;
    uu=uyz(:,m,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(:,m,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


wym=uzm;
wym_fft=uzmt_FFT;

disp(['Integral length scale wy: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,wym,'r-'); hold on
% plot(y_y,wym_fft,'r--')


%

%% wz
Nz=size(uyz,2);
Nlag=Nz-1;
uzm=0;
uzmt_FFT=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,3) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=uyz(m,:,k);%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
    
    uu_FFT = fft(squeeze(uyz(m,:,k)));
    uu_FFT(1)=0; % remove mean value. This is done inernally in function autocorr.
    uz_fft = ifft(abs(uu_FFT).^2);
    uz_fft=uz_fft/uz_fft(1);
    uzmt_FFT=uzmt_FFT+uz_fft;

end
end

uzm=uzm/cc/mm;
uzmt_FFT=uzmt_FFT/cc/mm;

iz=find(uzm<0,1);
Ly=trapz(y_y(1:iz),uzm(1:iz))*1000;


iz_fft=find(uzmt_FFT<0,1);
Lz_FFT=trapz(y_y(1:iz_fft),uzmt_FFT(1:iz_fft))*1000;


wzmm=uzm;
wzm_fft=uzmt_FFT;

disp(['Integral length scale wz: ',num2str(Ly),' mm',' and fft based: ' num2str(Lz_FFT) ' mm'])

% plot(y_y,wzmm,'b-'); hold on
% plot(y_y,wzm_fft,'b--')


%%

figure(1)
plot(y_y,uym_fft,'k')
hold on
plot(y_y,vym_fft,'b'); hold on

plot(y_y,wym_fft,'r')
plot(y_y,wzm_fft,'r--')
plot(y_y,vzm_fft,'b--')
plot(y_y,uzm_fft,'k--')

set(gca,'FontSize',16,'TickLabelInterpreter','latex')
ylabel('$ACF$','Interpreter','Latex','FontSize',18)
xlabel('$y$','Interpreter','Latex','FontSize',18)
legend('$uy$','$vy$','$wy$','$uz$','$vz$','$wz$','Interpreter','Latex')

figure(2)
plot(y_y,uym,'k'); hold on
plot(y_y,uym_fft,'k--')

plot(y_y,vym,'b'); hold on
plot(y_y,vym_fft,'b--'); hold on

plot(y_y,wym,'r'); hold on
plot(y_y,wym,'r--')

set(gca,'FontSize',16,'TickLabelInterpreter','latex')
ylabel('$ACF$','Interpreter','Latex','FontSize',18)
xlabel('$y$','Interpreter','Latex','FontSize',18)
legend('$uy$','$uy_{fft}$','$vy$','$vy_{fft}$','$wy$','$wy_{fft}$','Interpreter','Latex')

figure(3)

plot(y_y,uzmm,'k-'); hold on
plot(y_y,uzm_fft,'k--')

plot(y_y,vzmm,'b-')
plot(y_y,vzm_fft,'b--')

plot(y_y,wzmm,'r-'); 
plot(y_y,wzm_fft,'r--')

set(gca,'FontSize',16,'TickLabelInterpreter','latex')
ylabel('$ACF$','Interpreter','Latex','FontSize',18)
xlabel('$z$','Interpreter','Latex','FontSize',18)
legend('$uz$','$uz_{fft}$','$vz$','$vz_{fft}$','$wz$','$wz_{fft}$','Interpreter','Latex')

figure(4)
plot(y_y,uym,'k')
hold on
plot(y_y,vym,'b'); hold on

plot(y_y,wym,'r')
plot(y_y,uzmm,'k--')
plot(y_y,vzmm,'b--')
plot(y_y,wzmm,'r--')

set(gca,'FontSize',16,'TickLabelInterpreter','latex')
ylabel('$ACF$','Interpreter','Latex','FontSize',18)
xlabel('$y$','Interpreter','Latex','FontSize',18)
legend('$uy$','$vy$','$wy$','$uz$','$vz$','$wz$','Interpreter','Latex')
%%
% t_t=linspace(dt,5,Nt);

 uyz=ux;
Nt=size(uyz,3)-00000; % 20, 40, 60, 80, 90, 96, 98
t_final=2; % 40, 30, 20, 10, 5, 2, 1
t_t=linspace(dt,t_final,Nt);

%% ut
% Nt=size(uyz,3);
Nlag=Nt-1;
uzm=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,2) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=uyz(m,k,1:Nt);%-mean(uyz,'all'); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
end
end

uzm=uzm/cc/mm;
iz=find(uzm<0,1);
% y_y=linspace(0,0.15,Ny);
% t_t=linspace(0,1,Nt);
Lt=trapz(t_t(1:iz),uzm(1:iz))*Uinf*1000;
disp(['Integral length scale ut: ',num2str(Lt),' mm'])
%%
uyz=uy;

%% vt
Nt=size(uyz,3);
Nlag=Nt-1;
uzm=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,2) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=squeeze(uyz(m,k,:));%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
end
end

uzm=uzm/cc/mm;
iz=find(uzm<0,1);
% y_y=linspace(0,0.15,Ny);
% t_t=linspace(0,1,Nt);
Lt=trapz(t_t(1:iz),uzm(1:iz))*1000*Uinf;
disp(['Integral length scale vt: ',num2str(Lt),' mm'])
%%
uyz=uzz;

%% vt
Nt=size(uyz,3);
Nlag=Nt-1;
uzm=0;
mm=0;
cc = 0 ;
for k = 1:size(uyz,2) % loop over different time instances
    mm=mm+1;
    cc=0;
for m=1:size(uyz,1)
    cc = cc + 1 ;
    uu=squeeze(uyz(m,k,:));%-mean(u_1(k,:)); % In expeiments they may use the actual signal and not only the perturbation part
    uz = autocorr(uu,Nlag);
    uzm=uzm+uz;   
end
end

uzm=uzm/cc/mm;
iz=find(uzm<0,1);
% y_y=linspace(0,0.15,Ny);
% t_t=linspace(0,1,Nt);
Lt=trapz(t_t(1:iz),uzm(1:iz))*Uinf*1000;
disp(['Integral length scale wt: ',num2str(Lt),' mm'])
