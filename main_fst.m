%%
% close all
clear all
% clearvars -except seed
% clc
double precision;
format long;

addpath('./functions');
%% Parameters

% Calling parameters.m to get all parameters
parameters

disp(['ny= ',num2str(ny)])
disp(['nz= ',num2str(nz)])



%% Turbu ( -> spec -> sphere)

% calling turbu - generate wavenumbers and random amplitudes
turbu

%%
[Z,Y]=meshgrid(z_p,y_p);
Y=Y(:);
Z=Z(:);
X=Y*0+x_p(2);

kx=k_num_all(:,1);
ky=k_num_all(:,2);
kz=k_num_all(:,3);

% Phi_shift=bb1;%
Phi_shift=2*pi*rand(fst_modes,1);%

% Phi=X*kx' + Y*ky' + Z*kz';
% Phix=Phi+ ones(size(Y))*Phi_shift(:,1)'; % phase shifts different
% Phiy=Phi+ ones(size(Y))*Phi_shift(:,1)';
% Phiz=Phi+ ones(size(Y))*Phi_shift(:,1)';
Omega=Uinf*ones(size(Y))*kx';

Phi=X*kx' + Y*ky' + Z*kz' +ones(size(Y)).*Phi_shift' ;

%% compare with nek code
checks

%% Main time loop

U=struct;    % structures for saving velocity components

tic

for kt=1:N_tstep  % paralell loop used (for shorter tests use for)
    time=kt*dt;
    disp(['time=' num2str(time) ', timestep=' num2str(kt)])
    %     compute m sinus mode   
    
    s1 = sin( Phi + Omega*time );
%     s2 = sin( Phiy + Omega*time );
%     s3 = sin( Phiz + Omega*time );
    
    u1 =  s1*uamp';
    v1 =  s1*vamp';
    w1 =  s1*wamp';
 
        U(kt).u = reshape(u1,ny,nz,[]);
        U(kt).v = reshape(v1,ny,nz,[]);
        U(kt).w = reshape(w1,ny,nz,[]);
    
end

toc
% if ifsave
%  disp('Saving velocity')
%   save('Per-Mirr-6.mat','U','ym1_p','zm1_p','-v7.3')
%  save('fieldC24v.mat','Uv','ym1_p','zm1_p','-v7.3')
%  save('fieldC24w.mat','Uw','ym1_p','zm1_p','-v7.3')
%  disp('Script ended')
% end

% disp('Script ended')
%%
% Check field
ux=struct('U',{U.u});
ux=cell2mat(struct2cell(ux));

uy=struct('U',{U.v});
uy=cell2mat(struct2cell(uy));

uzz=struct('U',{U.w});
uzz=cell2mat(struct2cell(uzz));
%%
% plot fields
g=1;
figure(10)

pcolor(squeeze(zm1_p(1,1:end,1:end)),squeeze(ym1_p(1,1:end,1:end)),squeeze(ux(:,:,g)));
shading interp
colormap jet
axis equal
title('U')
colorbar
% pause

figure(11)

pcolor(squeeze(zm1_p(1,1:end,1:end)),squeeze(ym1_p(1,1:end,1:end)),squeeze(uy(:,:,g)));
shading interp
colormap jet
axis equal
colorbar
title('V')


figure(12)

pcolor(squeeze(zm1_p(1,1:end,1:end)),squeeze(ym1_p(1,1:end,1:end)),squeeze(uzz(:,:,g)));
shading interp
colormap jet
axis equal
colorbar
title('W')
hold on

figure(13)
urms=rms(ux-mean(ux,3),3);
vrms=rms(uy,3);
wrms=rms(uzz,3);
y=linspace(0,dly,ny);

  plot(y,mean(vrms,2)./mean(urms,2),'k-');hold on
plot(y,mean(wrms,2)./mean(urms,2),'k--')
xlabel('$y$','Interpreter','Latex')
ylabel('$v_{rms}/u_{rms}, w_{rms}/u_{rms}$','Interpreter','Latex')
set(gca, 'FontSize',18)
legend('$v_{rms}/u_{rms}$','$w_{rms}/u_{rms}$','Interpreter','Latex')

figure(14)
plot(y(1:end),3/2*(mean(urms,2).^2+mean(vrms,2).^2+mean(wrms,2).^2),'r-'); hold on
xlabel('$y$','Interpreter','Latex')
ylabel('$k$','Interpreter','Latex')
set(gca, 'FontSize',18)
legend('$k$','Interpreter','Latex')
ylim([0.4 0.6])
%%
Z_autoC

%%

um=mean(ux,3); um=mean(um,2);
vm=mean(uy,3); vm=mean(vm,2);
wm=mean(uzz,3); wm=mean(wm,2);

figure(56)
plot(y(1:end),um.^2,'r-'); hold on
plot(y(1:end),vm.^2,'g-');
plot(y(1:end),wm.^2,'b-');

