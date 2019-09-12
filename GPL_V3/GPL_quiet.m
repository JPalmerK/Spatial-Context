function [quiet_whiten, quiet_fft, quiet_base, noise_floor, blocked, baseline0] = GPL_quiet(sp,sp_whiten,parm)

low_off=parm.sum_bin_lo-parm.bin_lo+1;
high_off=parm.sum_bin_hi-parm.bin_lo+1;

u=bsxfun(@rdivide, sp_whiten, sum(sp_whiten.^2).^(1/2));
y=bsxfun(@rdivide, sp_whiten, sum(sp_whiten.^2, 2).^(1/2));

u=whiten_matrix(u')';
y=whiten_matrix(y);

bas=abs(u).^parm.xp1.*abs(y).^parm.xp2; 
baseline0=sum(bas(low_off:high_off,:).^2);

%%%%%%%%%%% Find quiet portions of slate

qs=sp_whiten;lks=0;quiet_fft=sp;
while lks < parm.nbin
    u=bsxfun(@rdivide, qs, sum(qs.^2).^(1/2));
    y=bsxfun(@rdivide, qs, sum(qs.^2, 2).^(1/2));
    bas=abs(u).^parm.xp1.*abs(y).^parm.xp2;
    baseline1=sum(bas.^2);

    ks=find(baseline1<=parm.noise_ceiling);
    qs=qs(:,ks);quiet_fft=quiet_fft(:,ks);

    lks=length(ks);
    xtnd=ceil(parm.nbin/lks);
    dex=zeros(1, xtnd*lks);
    for k=1:xtnd
        dex(lks*(k-1)+1:lks*k)=GPL_shuffle(lks,lks);
    end
    dex=dex(1:parm.nbin);
    qs=qs(:,dex);
    quiet_fft=quiet_fft(:,dex);
end

% Here is your quiet slate, 60 bins wide
quiet_fft=quiet_fft(:,1:60);

%%%%%%%%%%%%%%

%%%%%% First pass choice of quiet times

[b0,quiet_base]=whiten_vec(baseline0');
baseline0=b0'/quiet_base;
noise_floor=-min(baseline0);
ks=find(baseline0<=parm.noise_ceiling*noise_floor);
qs=sp_whiten(:,ks); 

% refine since when lots of calls, the small ones fall below noise level

tz2=size(qs, 2);
rep=ceil(parm.nbin/tz2);
qh=repmat(qs, 1, rep);
qs=qh(:,1:parm.nbin);

u=bsxfun(@rdivide, qs, sum(qs.^2).^(1/2));
y=bsxfun(@rdivide, qs, sum(qs.^2, 2).^(1/2));

u=whiten_matrix(u')';
y=whiten_matrix(y);

bas=abs(u).^parm.xp1.*abs(y).^parm.xp2;
baseline1=sum(bas(low_off:high_off,:).^2);

b1=whiten_vec(baseline1');

% scale with INITIAL mu_base
baseline1=b1'/quiet_base;
baseline1=baseline1(1:tz2);

ks_refine=find(baseline1>parm.noise_ceiling*noise_floor);
ks=ks(setdiff(1:tz2,ks_refine));

% be proactive - if below refined noise_ceiling at the start, disallow
blocked=zeros(parm.nbin,1);blocked(ks)=1;

qs=sp_whiten(:,ks);
tz2=size(qs,2);
rep=ceil(parm.nbin/tz2);
qh=repmat(qs, 1, rep);
quiet_whiten=qh(:,1:parm.nbin);
