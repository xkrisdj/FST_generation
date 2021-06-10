% random numbers on a sphere (radii - inside or 1 - on surface)
rvals = 2*rand(fst_modes,1)-1;
elevation = asin(rvals);

azimuth = 2*pi*rand(fst_modes,1);

radii = 1;%2*pi*(rand(fst_modes,1).^(1/3));

[bb1(:,1),bb1(:,2),bb1(:,3)] = sph2cart(azimuth,elevation,radii);

% get domain sizes
dlx=max(xm1_p(:))-min(xm1_p(:));
dly=max(ym1_p(:))-min(ym1_p(:));
dlz=max(zm1_p(:))-min(zm1_p(:));

% calling spec
spec; % get isotropically distributed wavenumbers in spheres 


% % using original rand function
% for k=1:ndim
%     bb(:,k) = ran2(seed)*2.0*pi;        % random phase shift 
%     bb1(:,k) = 2.0*ran2(seed)-1.0;      % random amplitude
% end

disp( 'FST - Random amplitude generated')
u_hat=zeros(fst_modes,ndim); 
u_hat_p=zeros(fst_modes,ndim);           
u_hat_pn=zeros(fst_modes,ndim);           

%     make sure that continuity is enforced by ensuring u_vec.k_vec=(0 0 0)
for i=1:fst_modes
    
  for j= 1:ndim
     u_hat(i,j)=bb1(i,j);
  end

  for j=1:ndim
     u_hat_p(i,j) = u_hat(i,j)...
         -(u_hat(i,1)*k_num_all(i,1)...
         + u_hat(i,2)*k_num_all(i,2)...
         + u_hat(i,3)*k_num_all(i,3))...
      * k_num_all(i,j)...
        / (k_num_all(i,1).^2 ...
        +  k_num_all(i,2).^2 ...
        +  k_num_all(i,3).^2);

  end
   
  for j=1:ndim
     u_hat_pn(i,j) = u_hat_p(i,j) ...
         / sqrt(u_hat_p(i,1).^2 ...
         + u_hat_p(i,2).^2 ...
         + u_hat_p(i,3).^2);
  end
end
disp('FST - Amplitudes projection done' )

%  Check energy in individual modes
ue=0.;
ve=0.;
we=0.;
uamp=zeros(1,fst_modes); 
vamp=zeros(1,fst_modes); 
wamp=zeros(1,fst_modes); 

for i=1:fst_modes
    shellno = shell(i);
    amp = shell_amp(shellno);
    uamp1 = u_hat_pn(i,1)*amp; %just for displaying
    vamp1 = u_hat_pn(i,2)*amp;
    wamp1 = u_hat_pn(i,3)*amp;

    ue =  ue + ((uamp1).^2)/2.; %just for displaying
    ve =  ve + ((vamp1).^2)/2.;
    we =  we + ((wamp1).^2)/2.;

    uamp(i) = u_hat_pn(i,1)*amp;
    vamp(i) = u_hat_pn(i,2)*amp;
    wamp(i) = u_hat_pn(i,3)*amp;
end
            
          
disp(['FST - Energy in u =' ' ' num2str(ue)]) 
disp(['FST - Energy in v =' ' ' num2str(ve) ])
disp(['FST - Energy in w =' ' ' num2str(we) ])
disp(['FST - Estimated tke =' ' ' num2str( (ue+ve+we)/2.)])
disp(['FST - Estimated ti =' ' ' num2str(sqrt((ue+ve+we)/3.))])
