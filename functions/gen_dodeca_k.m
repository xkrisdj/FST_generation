function [kx,ky,kz] = gen_dodeca_k(kk,Np,azimuth,elevation)
%GEN_DODECA_K Summary of this function goes here
%   Detailed explanation goes here

% 
% using original rand function
% rotx = ran2(seed)*2.*pi;
% roty = ran2(seed)*2.*pi;
% rotz = ran2(seed)*2.*pi;

% iy=1; iv=1;
% iseed1=-100;
% iseed2=-200;
% iseed3=-300;
% 
% rotx = randf_1(1,0,1,iseed1,iy,iv)*2.*pi;
% roty = randf_1(1,0,1,iseed2,iy,iv)*2.*pi;
% rotz = randf_1(1,0,1,iseed3,iy,iv)*2.*pi;
% 
rotx=rand(1,1)*2.*pi;
roty=rand(1,1)*2.*pi;
rotz=rand(1,1)*2.*pi;

% using random points on a sphere
% [rotx, roty, rotz] = sph2cart(azimuth,elevation,1);


% calling sphere
[kx,ky,kz]=sphere(Np,kk,rotx,roty,rotz);


end

