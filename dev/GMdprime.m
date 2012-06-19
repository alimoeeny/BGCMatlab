function d = GMdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
if size(G.mu,1) == 3
    %find three distances. one is allowed to be smaltt (splitting hash into
    %two. But one cluster must be distant from both of these. So take
    %lowest of top two = middle value
distance = mahal(G,G.mu);
    d(1) = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    d(2) = sqrt(2./((1./distance(1,3))+(1./distance(3,1))));
    d(3) = sqrt(2./((1./distance(2,3))+(1./distance(3,2))));
    ds = sort(d);
    d  = d(2);  
else
distance = mahal(G,G.mu);
    d = sqrt(2./((1./distance(2,1))+(1./distance(1,2))));
end
