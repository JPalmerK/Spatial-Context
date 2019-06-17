function idxN = nearest2(vector,refval)
% function to index the nearest value in a vector that is the same as the
% reference value.
%
% Where:
% vector is the vector of values you would like to index into
% refval is the referance value you want to know what is nearest to
% idxN is the index into the vector that is nearest to
%
% created by Dimitri W. Ponirakis
% Date: 12/17/08

idxN=[];

idxNh=find(vector>=refval);

if isempty(idxNh)==0

    idxNh=idxNh(1);
    
else
    
    idxNh=length(vector);
    
end


idxNl=find(vector<=refval);

if isempty(idxNl)==1
    
    idxNl=1;
    
end

idxNl=idxNl(end);

if idxNh==idxNl

    idxN=idxNh;

else

    if abs(vector(idxNh)-refval)<=abs(vector(idxNl)-refval)

        idxN=idxNh;

    else

        idxN=idxNl;

    end

end