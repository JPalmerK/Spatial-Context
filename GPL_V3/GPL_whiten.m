function sp_whiten = GPL_whiten(sp,parm)

[sp_whiten,mu]=whiten_matrix(sp,parm.white_x);
sp_whiten=abs(bsxfun(@rdivide, sp_whiten, mu));

% optional vertical whitening

if parm.whiten==2
    [sp_whiten,mu]=whiten_matrix(sp_whiten');
    sp_whiten=abs(sp_whiten'./(ones(parm.nfreq,1)*mu'));
end

cross=[0, 0.125, 0; 0.125, 0.500, 0.125; 0, 0.125, 0];
sp_whiten=conv2(sp_whiten,cross,'same');
