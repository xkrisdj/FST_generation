%
%     Philipp Schlatter's routine to generate wavenumbers on a sphere.
%
%---------------------------------------------------------------------- 

Nsmax=nshells;

il = fst_il;             
disp(['integral lenght scale =' ' ' num2str(il)])
Np = Npmax;
tke_scaled = (3.0/2.*fst_ti.^2);

%     Just initializing   
kxmax = 1.0E-20;
kxmin = 1.0E+20;

kymax = 1.0E-20;
kymin = 1.0E+20;

kzmax = 1.0E-20;
kzmin = 1.0E+20;

%     spectrum discretization
%     integrate the energy spectrum to determine total energy
Ndk = 5000;              % just a large no of points on the spectrum
dkint=(kend-kstart)/double(Ndk);
tke_tot1=(ek(kstart,il,1.)+ek(kend,il,1.)); 

for i=1:Ndk-1
    tke_tot1=tke_tot1+ek(kstart+i*dkint,il,1.); 
end
tke_tot1 = tke_tot1*dkint;

%     integrate the energy spectrum to determine total energy
%     with nshells      
dkint=(kend-kstart)/real(nshells-1);
tke_tot = 0.;

for i=1:nshells
    tke_tot=tke_tot+ek(kstart+(i-1)*dkint,il,1.); 
end

tke_tot = tke_tot*dkint;

disp(['FST - integrated energy in spectrum ' ' ' num2str(tke_tot1)])
disp(['FST - discretized on ' ' ' num2str(nshells) ' shells(used to scale) ' num2str(tke_tot)])


%     compute the coordinates using two dodecaeder

tke_tot1 = 0.;

for i=1:nshells

    k2=(kstart+(i-1)*(kend-kstart)/(nshells-1)).^2;
    kk(i)=sqrt(k2);
    dk(i) = (kend-kstart)/(nshells-1);
    q(i) = ek(sqrt(k2),il,tke_scaled/tke_tot);      % 1/tke_tot so the total 
                                                    % truncated energy = tke_scaled
    tke_shell(i) = q(i)*dk(i);
    tke_tot1 = tke_tot1 + tke_shell(i);

    % generating dodecaeder
    [cox(i,:),coy(i,:),coz(i,:)]=gen_dodeca_k(kk(i),Np,azimuth(i),elevation(i));

%     plot3(co1,co2,co3,'b.'); hold on

    % checking periodicity
    [cox(i,:),coy(i,:),coz(i,:)]=period(cox(i,:),coy(i,:),coz(i,:),Np,k2,dlx,dly,dlz,ifxp,ifyp,ifzp);
%     [cox(i,:),coy(i,:),coz(i,:)]=periodicity_chk(co1(i,:),co2(i,:),co3(i,:),Np,kk(i),dlx,dly,dlz,ifxp,ifyp,ifzp);
    


%     Get smallest and largest fst modes in x,y,z        
    ktmp1 = vlamax(cox(i,:),Np);
    if (ktmp1>kxmax) 
        kxmax = ktmp1;
    end
    ktmp2 = vlamin(cox(i,:),Np);
    if (ktmp2<kxmin) 
        kxmin = ktmp2;
    end
    ktmp3 = vlamax(coy(i,:),Np);
    if (ktmp3>kymax) 
        kymax = ktmp3;
    end
    ktmp4 = vlamin(coy(i,:),Np);
    if (ktmp4<kymin) 
        kymin = ktmp4;
    end
    ktmp5 = vlamax(coz(i,:),Np);
    if (ktmp5>kzmax) 
        kzmax = ktmp5;
    end
    ktmp6 = vlamin(coz(i,:),Np);
    if (ktmp6<kzmin) 
        kzmin = ktmp6;
    end
clear ktmp*
 end            % 1,nshells
 % mirroring
%  for i=1:nshells
% for j=Np+1:2*Np
% cox(i,j) = cox(i,j-Np);
% coy(i,j) = -coy(i,j-Np); %kdj: change from - to + to turn off mirroring
% coz(i,j) = -coz(i,j-Np); %kdj: change from - to + to turn off mirroring
% end
%  end
z1=0;
z2=0;
lu=zeros(1,nshells);
shell_modes=zeros(1,nshells);
l=0;

for i=1: nshells
    for j=1: Np
        %         If some modes need to be removed.
        %         Currently all modes except (0,0,0) are preserved.
        if (cox(i,j) == 0.0 && coy(i,j)==0.0 && coz(i,j)==0.0 )   
            z2=z2+1;
        else
            z1=z1+1;
            lu(i)=lu(i)+1;           % no of modes in each shell
            shell_modes(i)=shell_modes(i)+1;

            l=l+1;

            k_num_all(l,1) = cox(i,j);
            k_num_all(l,2) = coy(i,j);
            k_num_all(l,3) = coz(i,j);


            k_length =l;
            shell(l) = i;
        end          % if (.not.(0,0,0))
    end          % j=1,2*Np
end            % i=1,nshells

disp(['FST - (0,0,0) wavenumber removed'])
disp(['Saved ' num2str(z1) ' of ' num2str(z1+z2) ' fst modes.'])

tke_tot1 = 0.;
shell_amp=zeros(nshells,1);

for i=1:nshells
    shell_amp(i) = sqrt(2*tke_shell(i)*2./((shell_modes(i)))); % kdj check again
    
%      shell_amp(i) = sqrt(2*tke_shell(i)*dkint/((shell_modes(i)))); % kdj check again

    shell_energy = (shell_modes(i))*((shell_amp(i).^2.))/2.;
    tke_tot1 = tke_tot1 + shell_energy;
end

disp(['FST - Largest wavelength in x '  num2str(2.0*pi/kxmin)])
disp(['FST - Smallest wavelength in x ' num2str(2.0*pi/kxmax)])
disp(['FST - Largest wavelength in y '  num2str(2.0*pi/kymin)])
disp(['FST - Smallest wavelength in y ' num2str(2.0*pi/kymax)])
disp(['FST - Largest wavelength in z '  num2str(2.0*pi/kzmin)])
disp(['FST - Smallest wavelength in z ' num2str(2.0*pi/kzmax)])



