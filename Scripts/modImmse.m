function modErr = modImmse(mat1, mat2)

% Modified mean squared error for images. Excludes the area where both
% grids area equal
areaBothZero = intersect(find(mat1==0), find(mat2==0))


areaOneZero = unique([intersect(find(mat1==0), find(mat2 ~=0)),...
    intersect(find(mat1~=0), find(mat2 ==0))])
areaNeitherZero = intersect(find(mat1 ~=0), find(mat2 ~=0));

mat1(areaBothZero)=0/0;
mat2(areaBothZero)=0/0;




end