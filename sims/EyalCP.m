function resp = EyalCP(varargin)
%resp is mean resp for correct reject, false alarm, hits, misses
ntrials = 1000;
signals = 0.5;
noises = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noise',4)
        j = j+1;
        noises = varargin{j};
    elseif strncmpi(varargin{j},'signal',4)
        j = j+1;
        signals = varargin{j};
    elseif strncmpi(varargin{j},'ntrials',4)
        j = j+1;
        ntrials = varargin{j};
    end
j = j+1;
end

for k = 1:length(noises)
    noise = noises(k);
for j = 1:length(signals)
    signal = signals(j);
R(:,1) = randn(ntrials,1);
R(:,2) = randn(ntrials,1)+signal;
N = randn(size(R)).*noise;
D = R+N;
crit = signal/2;
x = -5:0.2:5;
dx=0.1;
hit = find(D(:,2) > crit);
miss = find(D(:,2) <= crit);
fa = find(D(:,1) > crit);
cr = find(D(:,1) <= crit);
a = hist(R(cr,1),x);
resp(1,j,k) = mean(R(cr,1));
bar(x,a);
hold on;
b = hist(R(fa,1),x,0.5);
resp(2,j,k) = mean(R(fa,1));
bar(x+dx,b,0.5,'r');
c = hist(R(hit,2),x);
resp(3,j,k) = mean(R(hit,2));

bar(x,-c,0.5);
d = hist(R(miss,2),x);
resp(4,j,k) = mean(R(miss,2));
bar(x+dx,-d,0.5,'r');
pc = 100*(length(hit)+length(cr))./length(R(:));
ratios(j,k) = (resp(3,j,k)-resp(2,j,k)+resp(4,j,k)-resp(1,j,k))./(2*signal);
fprintf('Hit-Fa %.2f, Miss-Reject %.2f %.1f correct Ratio %.2f\n',resp(3)-resp(2),resp(4)-resp(1),pc,ratios(j,k))
resp(5,j,k) = mean(R(:,1));
resp(6,j,k) = mean(R(:,2));
resp(7,j,k) = ratios(j,k);
resp(8,j,k) = (resp(3,j,k)-resp(2,j,k))./signal; %hits - FA. 
resp(9,j,k) = (resp(3,j,k)-resp(1,j,k))./signal; %hits - Correct Rejects. 
resp(10,j,k) = (resp(3,j,k)-resp(4,j,k))./signal; %hits - misses. 

end
end