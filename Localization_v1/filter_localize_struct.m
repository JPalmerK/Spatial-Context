function [localize_struct]=filter_localize_struct(localize_struct,array_struct,hyd,min_lsq_score);

parm=hyd(1).detection.parm;


for i=1:length(array_struct)
m(i)=array_struct(i).master;end

cc=[];
sc=[];rt=[];sr=[];
for k=1:length(m)
    
    localize_struct.hyd(m(k))
  % rt=localize_struct.hyd(m(k)).rtimes;
    
    n_dim=ndims(localize_struct.hyd(m(k)).coordinates);

    
 [s1,s2]=size(localize_struct.hyd(m(k)).score);
 size(localize_struct.hyd(m(k)).coordinates);
 if s2>0 
 if s1>1 && n_dim > 2
 localize_struct.hyd(m(k)).coordinates=squeeze(localize_struct.hyd(m(k)).coordinates(end-1,:,:));
 localize_struct.hyd(m(k)).score=localize_struct.hyd(m(k)).score(end-1,:);
 
 fg=localize_struct.hyd(m(k)).coordinates;
 req=localize_struct.hyd(m(k)).score;
 
 
 elseif s1 >0 && n_dim >2
 localize_struct.hyd(m(k)).coordinates=squeeze(localize_struct.hyd(m(k)).coordinates(1,:,:));
 localize_struct.hyd(m(k)).score=localize_struct.hyd(m(k)).score(1,:);
 
 fg=localize_struct.hyd(m(k)).coordinates;
 req=localize_struct.hyd(m(k)).score;
 
 else
 localize_struct.hyd(m(k)).coordinates= localize_struct.hyd(m(k)).coordinates(end-1,:)';
 localize_struct.hyd(m(k)).score=localize_struct.hyd(m(k)).score(end-1,:)';
 
 fg=localize_struct.hyd(m(k)).coordinates;
 req=localize_struct.hyd(m(k)).score;
 
 end
 %hom=localize_struct.hyd(m(k)).rtimes;
 xy=array_struct(k).ring;
 in=inpolygon(fg(1,:),fg(2,:),xy(:,1),xy(:,2));
 kk=find(in);
 
  kkk=find(req < min_lsq_score);
  
  kk=intersect(kk,kkk);
  
localize_struct.hyd(m(k)).delays=localize_struct.hyd(m(k)).delays(kk,:);
localize_struct.hyd(m(k)).cross_score=localize_struct.hyd(m(k)).cross_score(kk,:);
localize_struct.hyd(m(k)).coord_time=localize_struct.hyd(m(k)).coord_time(kk,:);
localize_struct.hyd(m(k)).rtimes=localize_struct.hyd(m(k)).rtimes(kk);
localize_struct.hyd(m(k)).dex=localize_struct.hyd(m(k)).dex(kk);
localize_struct.hyd(m(k)).coordinates=localize_struct.hyd(m(k)).coordinates(:,kk);
localize_struct.hyd(m(k)).score=localize_struct.hyd(m(k)).score(kk);
  
%   cc=localize_struct.hyd(m(k)).coordinates;
%   rt=localize_struct.hyd(m(k)).rtimes;
% 
% x=cc(2,:); y=cc(1,:); z=rt/1.e8;
% dd=[];
% 
% 
% for i=1:length(x)
%  dst=(x-x(i)).^2+(y-y(i)).^2+(z-z(i)).^2;
%  [k1,k2]=sort(dst);
%  if(k2 > 5)
%  dd(i)=dst(k2(6))^.5;
%  end
% end
% 
% 
% if(isempty(dd)==0)
% kkk=find(log(dd)<-5.5);
% 
% else
%     kkk=[];
% end
% 
% kk=kkk;
%  
%  
% localize_struct.hyd(m(k)).delays=localize_struct.hyd(m(k)).delays(kk,:);
% localize_struct.hyd(m(k)).cross_score=localize_struct.hyd(m(k)).cross_score(kk,:);
% localize_struct.hyd(m(k)).coord_time=localize_struct.hyd(m(k)).coord_time(kk,:);
% localize_struct.hyd(m(k)).rtimes=localize_struct.hyd(m(k)).rtimes(kk);
% localize_struct.hyd(m(k)).dex=localize_struct.hyd(m(k)).dex(kk);
% localize_struct.hyd(m(k)).coordinates=localize_struct.hyd(m(k)).coordinates(:,kk);
% localize_struct.hyd(m(k)).score=localize_struct.hyd(m(k)).score(kk);

end
end

