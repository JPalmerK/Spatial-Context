function [cs,zs] = effvp(c,z,zA,zB)
%---------------------------------------------------------------------
%
%    function [cs,zs] = effvp(c,z,zA,zB)
%
%   Produce sound velocity profile between  point A (zA) and point B (zB)
%
%   INPUT:
%
%   (c,z) = sound velocity profile (column vectors)
%   zA    = depth of point A
%   zB    = depth of point B
%
%   OUTPUT:
%
%    cs = sound velocity within the effective range (column vector)
%    zs = depth corresponding to cs
% 
%      by James Hu               7/12/95
%--------------------------------------------------------------------- 
%
if (zA > zB)
   tem = zA; zA = zB; zB=tem;
end
zm=max(z);
if (zm < zB),
    disp(['Maximum depth of SVP is ',num2str(zm)])
    zB=zm+20;
    disp(['This solution likely in error'])
    %error('Insufficient range (depth) for sound velocity profile')
end
%
%   restructure (cs,zs) working velocity profile
%
cA=interp1(z,c,zA);        % sound velocity at A
if zA>z(1),
   nA=max(find(zA>z));     % 
else
   nA=0;
end
nB = max(find(zB>z));
cB = interp1(z,c,zB);     % sound velocity at B
n = nB - nA +1;           % number of element
%
if zA==z(nA+1),
   cs=c(nA+1:nB);
   zs=z(nA+1:nB);
else
   cs=[cA;c(nA+1:nB)];
   zs=[zA;z(nA+1:nB)];
end
%
cs=[cs;cB]; zs=[zs;zB];
 
