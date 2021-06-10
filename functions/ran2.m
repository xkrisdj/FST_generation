function [ran21] = ran2(idum)
%
%     A simple portable random number generator
%
%     Requires 32-bit integer arithmetic
%     Taken from Numerical Recipes, William Press et al.
%     gives correlation free random numbers but does not have a very large
%     dynamic range, i.e only generates 714025 different numbers
%     for other use consult the above
%     Set idum negative for initialization
%
  
      m=714025;
      ia=1366;
      ic=150889;
      rm=1./m;
      iff=0;
      ir=zeros(97,1);

      if idum<0 | iff==0

%     Initialize

         iff=1;
         idum=mod(ic-idum,m);
         for j=1:97
            idum=mod(ia*idum+ic,m);
            ir(j)=idum;
         end
         idum=mod(ia*idum+ic,m);
         iy=idum;
      
      end
%     Generate random number

      j=1+floor(((97*iy)/m));
      iy=ir(j);
      ran21=iy*rm;
      idum=mod(ia*idum+ic,m);
      ir(j)=idum;
    
end
 
