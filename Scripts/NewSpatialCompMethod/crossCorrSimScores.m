function [dist, corrScore]= crossCorrSimScores( nextLklhdSpace, Lklhd_space_proj,...
         deltalat_space, deltalon_space)
     
     
     [I,J] = ind2sub(size(Lklhd_space_proj),find(Lklhd_space_proj>.05));
     if ~isempty(I)
         sxr =min(I):max(I);
         szc =min(J):max(J);
         
         
         % Use correlation to create the correlation map and find the translation
         %nimg = gridCall2-mean(mean(gridCall2));
         nimg = Lklhd_space_proj;
         nSec = nimg(sxr,szc);
         crr = xcorr2(nextLklhdSpace, nSec);
 
         [ssr,snd] = max(crr(:));
         [ij,ji] = ind2sub(size(crr),snd);
         
         
         %figure; imagesc(crr);
         
         normfac=sum(sum(nSec.^2));
         
         % The maximum of the cross-correlation corresponds to the estimated
         % location of the lower-right corner of the section. Use ind2sub to convert
         % the one-dimensional location of the maximum to two-dimensional coordinates.
         [ssr,snd] = max(crr(:));
         [ij,ji] = ind2sub(size(crr),snd);
         
         % Determin the spatial shift (x,y) between the maximum correlation lcoation
         % and the origional
         
         % lower right corner
         shiftr = (max(I)-(ij))*deltalat_space;
         shiftc = (max(J)-(ji))*deltalon_space;
         
         dist = sqrt(shiftr^2+shiftc^2);
         
         corrScore = ssr;
     else
         dist =Inf;
         corrScore =0;
         
 end