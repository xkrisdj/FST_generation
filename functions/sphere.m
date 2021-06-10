function [x,y,z] = sphere(Np,rad,rotx,roty,rotz)
  
%      this subroutine computes a set of Np point which are 
%      (more or less) uniformly distributed on a sphere with
%      radius R.
% 
%       input:    Np      Number of points
%                         (will be changed in case another number of
%                         points needs to be used)
%                 rad     Radius of sphere
%                 x,y,z   Arrays holding the coordinates (1..Np)
%                 rotx,roty,rotz  Rotation angles


N=1000;
new=0;
%--------------------------------------------------
if (Np>N)
    disp('number of points too big.')
    disp('Please change in sphere.f.')
return
end
if Np==20 && new==0
% disp('let us have a dodekaeder')
   
w=0.5*(sqrt(5.)+3.);
Np=20;
[xa1,ya1,za1]=asp(0.5*w,0.5*(w-1.),0.);
[xa2,ya2,za2]=asp(0.5*w,0.5*(w+1.),0.);
[xa3,ya3,za3]=asp(w-0.5,w-.5,.5);
[xa4,ya4,za4]=asp(w,.5*w,.5*(w-1.));
[xa5,ya5,za5]=asp(w-.5,0.5,.5);
[xa6,ya6,za6]=asp(.5*(w+1.),0.,.5*w);
[xa7,ya7,za7]=asp(.5*(w-1.),0.,.5*w);
[xa8,ya8,za8]=asp(.5,.5,.5);
[xa9,ya9,za9]=asp(0.,.5*w,(w-1.)*.5);
[xa10,ya10,za10]=asp(.5,w-.5,.5);
[xa11,ya11,za11]=asp(.5*(w-1.),w,.5*w);
[xa12,ya12,za12]=asp(.5*(w+1.),w,.5*w);
[xa13,ya13,za13]=asp(w-.5,w-.5,w-.5);
[xa14,ya14,za14]=asp(w,.5*w,.5*(w+1.));
[xa15,ya15,za15]=asp(w-.5,.5,w-.5);
[xa16,ya16,za16]=asp(.5*w,.5*(w-1.),w);
[xa17,ya17,za17]=asp(.5,.5,w-.5);
[xa18,ya18,za18]=asp(0.,.5*w,.5*(w+1.));
[xa19,ya19,za19]=asp(.5,w-.5,w-.5);
[xa20,ya20,za20]=asp(.5*w,.5*(w+1.),w);
    
xm=[xa1 xa2 xa3 xa4 xa5 xa6 xa7 xa8 xa9 xa10 xa11 xa12 xa13 xa14 xa15 xa16 xa17 xa18 xa19 xa20];
ym=[ya1 ya2 ya3 ya4 ya5 ya6 ya7 ya8 ya9 ya10 ya11 ya12 ya13 ya14 ya15 ya16 ya17 ya18 ya19 ya20];
zm=[za1 za2 za3 za4 za5 za6 za7 za8 za9 za10 za11 za12 za13 za14 za15 za16 za17 za18 za19 za20];

[xt1,yt1,zt1]=trans(xm,ym,zm,Np,-0.5*w,-0.5*w,-0.5*w);
[xs1,ys1,zs1]=scale1(xt1,yt1,zt1,Np,2./(sqrt(w.^2+1.)));

% next is not really used enywhere
Nl=30;
[l]=asl(1,1,2);
[l]=asl(2,2,3);
[l]=asl(3,3,4);
[l]=asl(4,4,5);
[l]=asl(5,5,6);
[l]=asl(6,6,7);
[l]=asl(7,7,8);
[l]=asl(8,8,9);
[l]=asl(9,9,10);
[l]=asl(10,10,11);
[l]=asl(11,11,12);
[l]=asl(12,12,13);
[l]=asl(13,13,14);
[l]=asl(14,14,15);
[l]=asl(15,15,16);
[l]=asl(16,16,17);
[l]=asl(17,17,18);
[l]=asl(18,18,19);
[l]=asl(19,19,20);
[l]=asl(20,3,12);
[l]=asl(21,19,11);
[l]=asl(22,13,20);
[l]=asl(23,16,20);
[l]=asl(24,7,17);
[l]=asl(25,6,15);
[l]=asl(26,1,5);
[l]=asl(27,1,8);
[l]=asl(28,2,10);
[l]=asl(29,4,14);
[l]=asl(30,9,18);

%     do the scaling to rad
% disp('doing scaling')
[xs,ys,zs]=scale1(xs1,ys1,zs1,Np,rad);

%     do the rotation
% disp('doing rotation')
[x,y,z]=rot3d(xs,ys,zs,Np,rotx,roty,rotz);


else
 %%   
 disp('Try new way')
k=0;
kk=0;

Nn=1;
Npn=0;
while Npn<=Np
Npn=0;
dtheta = 2*pi / real(Nn);
for j=0: (Nn)/2
    theta=j*dtheta;
    if (sin(theta) == 0.0)
        dphi = 99999999.;
    else
        dphi = dtheta / sin(theta);
    end
    Nphir = max(2*pi / dphi ,1.);
    Nphi = Nphir+0.5;
    dphi = 2*pi / real(Nphi);
    Npn=Npn+Nphi;
end

    Npn1=Npn;
    Nn=Nn+1;
end
Npn=Npn1;
Nn=Nn-1;

dtheta = 2*pi / real(Nn);

for j=0: (Nn)/2
    theta=j*dtheta;
    if (sin(theta)==0.0)
        dphi = 99999999.;
    else
        dphi = dtheta / sin(theta);
    end
    Nphir = max(2*pi / dphi ,1.);
    Nphi = Nphir+0.5;
    dphi = 2*pi / real(Nphi);
    for i=1: (Nphi)
        phi=i*dphi;
        k=k+1;
        xm(k) = cos(phi)*sin(theta);
        ym(k) = sin(phi)*sin(theta);
        zm(k) = cos ( theta );
        kk=kk+1;
        l=asl(kk,k,k+1);
    end
end
Np=k;
Nl=kk-1;

%     do the scaling to rad
% disp('doing scaling')
[xs,ys,zs]=scale1(xm,ym,zm,Np,rad);

%     do the rotation
% disp('doing rotation')
[x,y,z]=rot3d(xs,ys,zs,Np,rotx,roty,rotz);

end
%%



% %       check...
%       for i=1: Np
%         d=x(i).^2.+y(i).^2.+z(i).^2.;
%         disp([num2str(d) ' ' num2str(d.^2.)])
%       end 


end

