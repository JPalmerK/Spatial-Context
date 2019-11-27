function drawSimMat(obj, titlestr)
% Draw the simualtion class and make it look nice

figure;

imAlpha=ones(size(obj.Sim_mat));
imAlpha(isnan(obj.Sim_mat))=0;
imagesc(obj.Sim_mat,'AlphaData',imAlpha);
set(gca,'color',[1 1 1]);


% 
% 
% 
% [nr,nc] = size(obj.Sim_mat);
% h =pcolor(flipud((obj.Sim_mat)));
% set(h, 'EdgeColor', 'none');
% ax = gca;
% ax.YTickLabel = flipud(ax.YTickLabel)
% colormap(jet);
% colorbar
% 
% xlabel('Call ID')
% ylabel('Call ID')
% if nargin==2
%     title(obj.titleStr)
% end


end