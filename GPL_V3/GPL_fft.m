function [sp] = GPL_fft(data,parm);

dt=parm.skip/parm.sample_freq;
cutoff_short=ceil(parm.min_call/dt)-2;
cutoff_long=ceil(parm.max_call/dt)-2;
win=hamming(parm.fftl);
sp=zeros(parm.nfreq,parm.nbin);

x=data;[x1,x2]=size(x); 
if x2>x1
  x=x';end

% error condition causing out of bounds fft operation
% parm.nbin needs to match such that fft(x(start:finish)) meets the
% crierion mod(length(x), (finish-start-1) == 0.

% % disp('==================');

for j=1:parm.nbin
    
    try
        start=(j-1)*parm.skip+1;
        finish=start+parm.fftl-1;
        xl = length(x(start:finish));
        wl = length(win);
        q=fft(x(start:finish).*win);  % errors here.              
        sp(:,j)=abs(q(parm.bin_lo:parm.bin_hi));
    catch
        % something bad happened, likely a incomplete sound series, return.
%         disp('Error : Break');
        sp = [];
        return;
    end
    
end
