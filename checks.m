% check spectra
figure(2)
loglog(sqrt(kx.^2+ky.^2+kz.^2),sqrt(uamp.^2+vamp.^2+wamp.^2),'r-','LineWidth',1.7); hold on
% plot von Karman spectra
% loglog
% plot spectra from nek code
% fst=importdata('../fst_spectrum.dat');
% k1=fst.data(:,2);
% k2=fst.data(:,3);
% k3=fst.data(:,4);
% 
% E1=fst.data(:,5);
% E2=fst.data(:,6);
% E3=fst.data(:,7);

% k_nek=sqrt(k1.^2+k2.^2+k3.^2);
% E_nek=sqrt(E1.^2+E2.^2+E3.^2);
% loglog(k_nek,E_nek,'b--','LineWidth',2)

hold off
legend('$E_{matlab}$','$E_{nek}$','Interpreter','Latex')

xlabel('$k$','Interpreter','Latex');
ylabel('$E$','Interpreter','Latex');
set(gca,'FontSize',18,'TickLabelInterpreter','latex')

% check how wavenumbers are distributed
figure(1)
plot3(kx,ky,kz,'r.'); hold on
% plot3(k1,k2,k3,'b.')
xlabel('$k_x$','Interpreter','Latex')
ylabel('$k_y$','Interpreter','Latex')
zlabel('$k_z$','Interpreter','Latex')
axis equal
% % check random phase
% load('bb_nek.txt');
% figure(3)
% plot(bb(:,1),'r.'); hold on
% plot(bb(:,2),'b.'); hold on
% plot(bb(:,3),'g.'); hold on
% 
% bb_nek=reshape(bb_nek,[],6);
% plot(bb_nek(:,1),'ro'); hold on
% plot(bb_nek(:,3),'bo')
% plot(bb_nek(:,5),'go')