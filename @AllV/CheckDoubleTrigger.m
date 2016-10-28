function newTrigger = CheckDoubleTrigger(DATA, AllVoltages, varargin);
%AllC.CheckDoubleTrigger detects when a spike waveform is causing double
%triggers

id = varargin{1};
sid = find(diff(id) < 10);
e = DATA.energy(1,sid);
if mean(e) > prctile(DATA.energy(1,:),50);
    DATA = AllV.AddErr(DATA,'Double Triggers are Real Spikes\n');
end
