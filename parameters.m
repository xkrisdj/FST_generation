%% Parameters

nshells = 180;        % No of spherical shells
Npmax = 20;       	 % Shell discretization points

fst_ti = 0.051*6.0;  % turbulence intensity
fst_il = 2.5E-02;    % integral length scale
	
kstart = 42.00;      % smallest wavenumber

kend = 1660.;        % largest wavenumber 	

ndim=3; % number of dimensions
Uinf=6; % translation velocity


dt=5e-4; % timestep
t_final=2.; % final time in sec
N_tstep=t_final/dt; % number of iterations


% Periodicity
ifxp = 0; %if periodic in x
ifyp = 0; %if periodic in y
ifzp = 0; %if periodic in z



%% Initialize

fst_modes=nshells*Npmax; % No of freestream modes
shell_modes=zeros(nshells,1); % Modes saved per shell    

bb = zeros(fst_modes,3); % random phase
bb1 = zeros(fst_modes,3); % random amplitude

shell=zeros(fst_modes,1);
shell_amp=zeros(nshells,1);

k_num_all=zeros(fst_modes,3); % all wavenumbers
kk=zeros(nshells,1);

u_hat=zeros(fst_modes,3);
u_hat_p=zeros(fst_modes,3);
u_hat_pn=zeros(fst_modes,3);

%% Create mesh xm1_p,ym1_p,zm1_p
xmax=0.1; xmin=-0.1;
ymax=xmax; ymin=xmin;

ny=128; nz=128;

x_p=linspace(-0.183,1.165,5);                  % creating mesh
y_p=linspace(xmin,xmax,ny); 
z_p=linspace(ymin,ymax,nz); 

[xm1_p, ym1_p, zm1_p]=ndgrid(x_p,y_p,z_p);

