function [spc, mu]=whiten_matrix(sp,fac)

if nargin==1
    fac=1;
end
mu=base3x(sp);
if fac == 0
    spc=sp;
else
    spc=bsxfun(@minus, sp, mu*fac);
end
