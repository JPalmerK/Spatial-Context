function mu = base3x(baseline)

if ~isvector(baseline)
    bs=sort(baseline, 2);
    [rows, cols]=size(baseline);
    bs1=floor(cols/2);
    [~, b2]=min(bs(:, bs1+1:2*bs1)-bs(:, 1:bs1), [], 2);
    ix = repmat((0:bs1)*rows, rows, 1);
    offset = (0:rows-1)' + (b2-1)*rows + 1;
    ix = ix + offset;
    qs = bs(ix);
    mu = qs*ones(bs1+1, 1)/(bs1+1);
else
    bs = sort(baseline);
    bs1=floor(length(bs)/2);
    [~, b2]=min(bs(bs1+1:2*bs1)-bs(1:bs1));
    mu=sum(bs(b2:b2+bs1))/(bs1+1);
end
