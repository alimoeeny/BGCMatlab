function Expt = CheckEyeData(Expt, varargin)
%Expt = CheckEyeData(Expt) makes sure lengths of EyeData
%in Trials match. 

zeropad = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',5)
        zeropad = 1;
    end
    j = j+1;
end


for j = 1:length(Expt.Trials)
    emlen(j) = length(Expt.Trials(j).EyeData);
    if emlen(j) > 0
        if isnan(Expt.Trials(j).emstarti)
            prepts(j) = 0;
        else
            prepts(j) = Expt.Trials(j).emstarti;
        end
        postpts(j) = emlen(j)-prepts(j);
    else
        prepts(j) = 0;
        postpts(j) = 0;
    end
end

pre = max(prepts);
post = max(postpts);
for j = 1:length(Expt.Trials)
    sz = size(Expt.Trials(j).EyeData,2); %might have 4,6 or 7 channels
    prepad = ones(pre-prepts(j),sz).*NaN;
    postpad = ones(post-postpts(j),sz).*NaN;
    E = cat(1, prepad, Expt.Trials(j).EyeData, postpad);
    if zeropad
        E(isnan(E)) = 0;
    end
    Expt.Trials(j).EyeData = E;
    if isfield(Expt.Trials,'pupil')
        E = cat(1, prepad(:,1), Expt.Trials(j).pupil(:), postpad(:,1));
        Expt.Trials(j).pupil = E;
    end
end
emlen = post+pre;
Expt.Header.emlen = post+pre;
Expt.Header.emtimes = (([1:emlen]-pre) .* 1./Expt.Header.CRsamplerate);

min(emlen);