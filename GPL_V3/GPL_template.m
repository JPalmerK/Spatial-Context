function [GPL_struct] = GPL_template(GPL_struct,sp,sp_whiten,start,finish,x,parm);

[m1,m2]=size(start);

for k=1:m1
    
base=sp_whiten(:,start(k):finish(k)); 
[cm,cm_max,cm_max2]=GPL_contour(base,parm.cut);

if(parm.cm_on == 1)
GPL_struct = GPL_sparse(cm,'cm',k,GPL_struct);
end
if(parm.cm_max_on==1)
GPL_struct = GPL_sparse(cm_max,'cm_max',k,GPL_struct);
end
if(parm.cm_max2_on==1)
GPL_struct = GPL_sparse(cm_max2,'cm_max2',k,GPL_struct);
end

% reconstruct with:
%  cm = GPL_full('cm',k,GPL_struct);

if parm.waveform==1
%  xv=x(1+(start(k)-1)*skip:(finish(k)-1)*skip);
%  [vec]=ww365(xv,cm_max,parm,sp(:,start(k):finish(k)));
%  GPL_struct(k).waveform=vec;end

%Modified on 2/23/19 by T.Helble to save waveform attributes instead of
%GPL spectral waveform
xv=x(1+(start(k)-1)*parm.skip:(finish(k)-1)*parm.skip);
xvf=bandpass(xv,[1250 1600],parm.sample_freq);

GPL_struct(k).waveform_rln_dB = 20*log10(sqrt(mean(xvf.^2)));
%GPL_struct(k).waveform_sel_dB = 20*log10(sqrt(mean(xvf.^2))) +
%10*log10(length(xvf)/parm.sample_freq); calculate this later when noise is
%removed
GPL_struct(k).waveform_pp_rln_dB = 20*log10(max(xvf) - min(xvf)); 
end

end