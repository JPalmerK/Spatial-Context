function [err,delay_p]=nelder_mead_2D(vec,array,delay,dep,svp,zz);
% freeze depth at 5 m
%E = [6378137.0 0.081819190842613];
whale=vec;

[sz1,sz2]=size(array);
n=length(array)-1;
for j=1:sz1
    % Matlab geodesic on WGS84 ellipsoid 
    hor=vdist(whale(1),whale(2),array(j,1),array(j,2));
    [s1,s2,speed] = aute(svp,zz,dep(j),hor,5);
    t(j)=((dep(j)-5)^2+hor^2)^.5/speed;end;
delay_p=t(2:end)-t(1)*ones(1,n);
v=delay_p-delay;
tt0=inv(.4*ones(n,n)+0.6*eye(n));
err=1-exp(-v*tt0*v');

