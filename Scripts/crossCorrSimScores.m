function [dist, corrScore]= crossCorrSimScores(obj, nextLklhdSpace, Lklhd_space_proj,...
         deltalat_space, deltalon_space)
     % X and Y values of ambiguity surface 1 where values >0
     
        [I,J] = ind2sub(size(Lklhd_space_proj),find(Lklhd_space_proj>0));

         sxr =min(I):max(I);
         szc =min(J):max(J);
         
         
         % Use correlation to create the correlation map and find the translation
         %nimg = gridCall2-mean(mean(gridCall2));
         nimg = Lklhd_space_proj;
         nSec = nimg(sxr,szc);
         
         crr = xcorr2(nextLklhdSpace, Lklhd_space_proj);

         
         [ssr,snd] = max(crr(:));
         [ij,ji] = ind2sub(size(crr),snd);
         
         
         %figure; imagesc(crr);
         
         
         % The maximum of the cross-correlation corresponds to the estimated
         % location of the lower-right corner of the section. Use ind2sub to convert
         % the one-dimensional location of the maximum to two-dimensional coordinates.
         [ssr,snd] = max(crr(:));
         [ij,ji] = ind2sub(size(crr),snd);
         
         % Determin the spatial shift (x,y) between the maximum correlation lcoation
         % and the origional
         
         % lower right corner
         shiftr = max([0, (max(I)-(ij+1))*deltalat_space]);
         shiftc = max([0,(max(J)-(ji+1))*deltalon_space]);
         
         dist = sqrt(shiftr^2+shiftc^2);
         
         corrScore = ssr;

 end