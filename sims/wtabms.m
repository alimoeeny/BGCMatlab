function wtabms(varargin)
%
%implement Winner take all with an output nonlinearity, and see how
%this models responses to Cuts

method = 'opnl';

x = -60:60;
sigma = 10;
G = Gauss(sigma,x);
plot(x,G);
gamma  = 2;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'gamma',5)
        j = j+1;
        gamma = varargin{j};
    elseif strncmpi(varargin{j},'divis',5)
        method = 'divis';
    end
    j = j+1;
end


if strcmp(method,'opnl')
%G is a gaussian response to SF
%G^gamma is response to these componenets after output NL
R(1,:) = (G).^gamma;
%sum(G)^gamma is response to the sum - a broadband stimulua
% (sum(G) -G).^gamma is the response to broadband with cuts
R(2,:) = (sum(G) -G).^gamma;
%invert and normalize these to compare with components alone
R(3,:) = -(R(2,:)-max(R(2,:)));
R(3,:) = R(3,:) .* max(R(1,:))/max(R(3,:));
plot(x,R);
legend('Response to components','Response to Cuts','Cuts compared');

elseif strcmp(method,'divis')
    for j = 1:length(gamma)
        g = gamma(j);
        R(1,:) = (G);
        %sum(G)^gamma is response to the sum - a broadband stimulua
        % (sum(G) -G).^gamma is the response to broadband with cuts
        R(2,:) = (sum(G.^g) -G.^g).^(1/g);
        %invert and normalize these to compare with components alone
        R(3,:) = -(R(2,:)-max(R(2,:)));
        R(3,:) = R(3,:) .* max(R(1,:))/max(R(3,:));
        R(2,:) = R(2,:) .* max(R(1,:))/max(R(2,:));
        plot(x,R);
        legend('Response to components','Response to Cuts','Cuts compared');
        hold on;
        allr(j,:) = R(2,:) ./max(R(2,:));
        
        for k = 1:length(G)
            sx(k) = (k * 0.5)-length(G)/4;
            sG = Gauss([ sx(k) 0.5 * sigma],x);
            sG = sG./max(sG);  %1 at peak
            sR = sG .* G; %response to components
            xR = (1-sG) .*G;  %response to cut components
            R(4,k) = sum(sR.^g).^(1/g);
            R(5,k) = sum(xR.^g).^(1/g);
        end
        R(4,:) = R(4,:) .* max(R(1,:))/max(R(4,:));
        R(5,:) = R(5,:) .* max(R(1,:))/max(R(5,:));
        plot(sx,R(5,:));
        R(5,:) = (max(R(5,:)) - R(5,:));
        R(5,:) = R(5,:) .* max(R(1,:))/max(R(5,:));
        plot(sx,R(4,:));
        plot(sx,R(5,:));
    end
    if length(gamma) > 1
        hold off;
        plot(allr');
    end
end