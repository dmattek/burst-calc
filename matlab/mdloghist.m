function [x y edges] = mdloghist(data, factor)

% histogram with logarithmic binning
% the size of the next edge is 'factor' bigger than the previous one
% Input:
% data - vector with data for binning
% factor - factor of bin increase
%
% Output:
% x - positions of the bins (at half of the width)
% y - normalized counts
% edges - position of the edges
%
% Example usage:
% inter = load('output1.prel');
% [x y edges] = mdloghist(inter, 1.5);
% loglog(x,y);
%

% the minimum number in the data vector
xmin=min(data);

% automatically computes the number of bins: nbins
% formula: xmin * factor^nbins = max
nbins = ceil(log(max(data)/xmin) / log(factor));

% computing edges
edges = zeros(nbins, 1);
edges(1) = xmin;
for i = 2:nbins
    edges(i) = edges(i-1)*factor;
end

% arguments for histogram
x = edges(1:(nbins-1))+(diff(edges)./2);

% binning data
dp = histc(data, edges);
ry = dp(1:(nbins-1));
dedges = diff(edges);

% normalization
y = ry./(sum(ry).*dedges);

% print the number of bins
nbins
