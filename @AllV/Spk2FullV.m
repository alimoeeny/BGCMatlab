function V = Spk2FullV(DATA, t, varargin)
%Spk2FullV(DATA, t, varargin) Build a dummy FullV from Spks
%

A = getappdata(DATA.toplevel,'AllVoltages');
id = find(DATA.t > t(1) & DATA.t < t(2));

for j = 1:length(id)
    s = round((DATA.t(id(j))-t(1)) * DATA.samplerate);
    v(s:s+size(A,2)-1) = A(DATA.probe(1),:,id(j));
end
V.t = t(1):1./DATA.samplerate:t(2);
v(end+1:length(V.t)) = 0;
V.V(DATA.probe(1),:) = v;
