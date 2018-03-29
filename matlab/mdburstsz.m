function [tn cn] = mdburstsz(data, maxsz, fac)

% computes burst size 'Phi(tarr)' as a function of threshold 'tarr'
% formula: Phi (tarr) = total no. of arrivals / number of arrivals larger than tarr
%
% data - input intervals
% maxsz - the total number of thresholds
% fac - factor to scale the interval range (type as 1/x to have x equal maximum threshold size)
%
% Example usage:
% inter = load('output1.prel');
% [t c] = mdburstsz(inter, 100, 1/2);
% loglog(t,c);
%


N = length(data)
maxth = 1/(fac*maxsz);

tn=zeros(maxsz, 1);
cn=zeros(maxsz, 1);

for i=1:maxsz
    nsmall=length(find(data>i*maxth));
    bsz=N/nsmall; 
    cn(i,1)=bsz;
    tn(i,1)=i/(fac*maxsz); 
end
