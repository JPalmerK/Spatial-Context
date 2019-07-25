function [corrScore]= crossCorrSimScoresFilt(AmbigNext, AmbProj, ...
    filtLat, filtLon)
% X and Y values of ambiguity surface 1 where values >0

[I,J] = ind2sub(size(AmbProj),find(AmbProj>.5));
sxr =min(I):max(I);
szc =min(J):max(J);

sxrb =min(I)-filtLat:max(I)+filtLat;
szcb =min(J)-filtLon:max(J)+filtLon;

% Use correlation to create the correlation map and find the translation
%nimg = gridCall2-mean(mean(gridCall2));
nimg = AmbProj;
nSec = nimg(sxr,szc);


sxrb =sxrb(logical(sxrb>0));
szcb =szcb(logical(szcb>0));


sxrb =sxrb(logical(sxrb<size(AmbigNext,1)));
szcb =szcb(logical(szcb<size(AmbigNext,2)));

imageB = AmbigNext(sxrb, szcb);
%imagesc(imageB)


crr = xcorr2(imageB, nSec);


[ssr,snd] = max(crr(:));
[ij,ji] = ind2sub(size(crr),snd);


%figure; imagesc(crr);

normfac=sum(sum(nSec.^2));

% Normalize the correlation
%crr = crr/normfac;

% The maximum of the cross-correlation corresponds to the estimated
% location of the lower-right corner of the section. Use ind2sub to convert
% the one-dimensional location of the maximum to two-dimensional coordinates.
ssr = (max(crr(:)));


corrScore = ssr;


end