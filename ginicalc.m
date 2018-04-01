function g = ginicalc(Wrow)

[~,ndyn] = size(Wrow);
Wsort = sort(Wrow,'ascend');
Wsum = sum(Wsort);
bigW = repmat(Wsort,ndyn,1);
g = (ndyn-1)/ndyn - 2*sum(sum(tril(ones(ndyn),-1).*bigW))/(ndyn*Wsum);