function [ve2,ve3,ve4] = aute(svp,zz,depth,distance,whale);
% -----------------------------------------------------------
%  function ve = aute(svp,zz,depth,distance);
%
%
%  Calculate effective mean sound velocity according to the
%  current AUTEC algorithm.
% 
% 
%
% INPUT:
%  svp : stratified sound velocity profile (column vector)
%  zz : depth corresponding to svp (column vector)
%  depth : depth for hydrophone
%  distance : range from whale
%
% OUTPUT:
%
%  ve2 : effective sound velocity (second-order series expansion)
%
%                James Hu, 7/29/95
%
% Ref: Vass, A.E. (1964) "Refraction of the direct monotonic 
%      sound ray", TM No. 322, p.3.
%---------------------------------------------------------------
%

if nargin==4
whale=5;end

[v,z] = EFFVP(svp,zz,whale,depth);  % effective range for SVP
vd =depth-whale;        % vertical distance
vbar=trapz(z,v) / vd;   % arthimetic mean sound velocity
ctb = distance/vd;
C2 = 1 - 1/2*ctb^2;
intgr2 = ((v-vbar)/vbar).^2;
J2 = trapz(z,intgr2) / vd;
fac1 = 1-C2*J2;
ve2=vbar*fac1;

if (nargout >1)
  C3 = -(C2 + 0.58*ctb^4);
  C4 = (ctb^4 /8) *(3-5*ctb^2) + C2; 
  intgr3 = ((v-vbar)/vbar).^3;
  intgr4 = ((v-vbar)/vbar).^4;
  J3 = trapz(z,intgr3) / vd;
  J4 = trapz(z,intgr4) / vd;
  fac2=fac1-C3*J3;
  ve3=vbar*fac2;
  fac3=fac2-C4*J4;
  ve4=vbar*fac3;
end
