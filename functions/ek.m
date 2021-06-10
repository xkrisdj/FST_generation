      function [ek]= ek(k,L,q)
%
%     gives the energy spectrum according to von Karman
%     input  k    wavenumber
%            L    integral length scale
%            q    total turbulent kinetic energy
%
%     ek statisfies the following conditions
%         1.)  integral of ek over k=0..infty is q (tot.kin.energy)
%         2.)  maximum of ek is at k=1.8 / L (see Tennekes p.272)
%
     
      ek = 2./3.*q*1.606 * (k*L).^4. * L /...
          (1.350+(k*L).^2.).^(17./6.);
      end
