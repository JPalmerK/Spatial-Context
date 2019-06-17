function [err,delay_p]=postfix2gls(vec,vec0,array,delay,dep,svp,zz,rho);

E = [6378137.0 0.081819190842613];
whale(1:2)=vec0(1:2)+vec(1:2)/1.e5;
whale(3)=vec0(3)+vec(3);

for j=1:5
    % Matlab geodesic on WGS84 ellipsoid 
    hor=distance(whale(1),whale(2),array(j,1),array(j,2),E);
    [s1,s2,speed] = aute(svp,zz,dep(j),hor,whale(3));
    t(j)=((dep(j)-whale(3))^2+hor^2)^.5/speed;end;
delay_p=t([1,2,4,5])-t(3)*ones(1,4);
v=delay_p-delay;
if sum(size(rho))==2
tt0=inv(rho*ones(4,4)+(1-rho)*eye(4));
else 
    tt0=inv(rho);end

err=1-exp(-v*tt0*v');

