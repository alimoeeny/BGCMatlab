function ddm(slopes, crit)
%ddm(slopes, crit)  Quick Drift diffusion model
%  ddm([0.05:0.05:0.4], 5) is default
if nargin == 0
    slopes = 0.05:0.05:0.4;
end
if nargin < 2
    crit = 5;
end
ntrials = 10000;
nts = 1000;

for j = 1:length(slopes);
clear rt;
clear score;
    rnd = randn(ntrials, nts);
    slope = [1:nts] .* slopes(j);
    signal = repmat(slope,ntrials,1);
    crnd = cumsum(rnd,2);
    resp = crnd+signal;
    for t = 1:ntrials
        a = find(abs(resp(t,:)) > crit,1);
        if isempty(a)
            score(t) = resp(t,end) > 0;
            rt(t) = nts+1;
        else
            score(t) = resp(t,a) > 0;
            rt(t) = a;
        end
    end
    pcorrect(j) = mean(score);
    rts(j) = mean(rt);
end
subplot(2,1,1);
plot(slopes,pcorrect);
subplot(2,1,2);
plot(slopes,rts);
%imagesc(resp);