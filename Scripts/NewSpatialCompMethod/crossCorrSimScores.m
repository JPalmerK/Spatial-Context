function [dist, corrScore]= crossCorrSimScores(gridCall1,gridCall2, deltalat_space, deltalon_space)
% X and Y values of ambiguity surface 1 where values >0

[I,J] = ind2sub(size(gridCall1),find(gridCall1>.1));
White = max(max(gridCall1));
sxr =min(I):max(I);
szc =min(J):max(J);

%normfac = max(length(find(gridCall1>.1)),length(find(gridCall2>.1)));


Sect = gridCall1(sxr, szc);
kimg = gridCall1;
kimg(sxr, szc) = White;


kumg = White*ones(size(gridCall1));
kumg(sxr, szc) = Sect;


% Use correlation to create the correlation map and find the translation
%nimg = gridCall2-mean(mean(gridCall2));
nimg = gridCall2;
nimg1 = (nimg-mean(nimg(:)))/((std(nimg(:)) * prod(size(nimg))));

nSec = gridCall1(sxr,szc);
nSec1 = (nSec-mean(nSec(:)))/((std(nSec(:)) * prod(size(nSec))));
crr = xcorr2(nSec, nimg);
%figure; imagesc(crr);



% The maximum of the cross-correlation corresponds to the estimated
% location of the lower-right corner of the section. Use ind2sub to convert
% the one-dimensional location of the maximum to two-dimensional coordinates.
%R = corr2(gridCall1,gridCall2);
[ssr,snd] = max(crr(:));
[ij,ji] = ind2sub(size(crr),snd);

% Determin the spatial shift (x,y) between the maximum correlation lcoation
% and the origional 

% lower right corner
shiftr = max([0, (max(I)-ij)*deltalat_space]);
shiftc = max([ 0, (max(J)-ji)*deltalon_space]);

dist = sqrt(shiftr^2+shiftc^2);

corrScore = snd;


end