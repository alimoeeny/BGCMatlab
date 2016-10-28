 function PlotMenu(a, b, type, varargin)      DATA = GetDataFromFig(a);     onoff = {'off' 'on'};     if sum(strcmp(type,{'xyseq' 'xcorrprobes'}))         DATA.plot.(type) = ~DATA.plot.(type);         set(a,'Checked',onoff{DATA.plot.(type)+1});         if DATA.plot.(type)             if strcmp(type,'xyseq')                 AllV.PlotXYSequence(DATA,DATA.cluster);             elseif strcmp(type,'xcorrprobes')                 AllV.CalcXcorr(DATA,[],'probes');             end         end     elseif sum(strcmp(type,'xcorrNone'))         DATA.plot.xcorrprobes = 0;         DATA.plot.xcorrtype = 'None';         AllV.SetGUI(DATA);     elseif sum(strcmp(type,'meanmenu'))         if sum(strcmp(varargin{1},{'sidebyside' 'image+lines' 'dprimeimage' 'allclusters'}))             DATA.plot.meantype = varargin{1};             DATA.plot.comparemean = 0;             AllV.PlotMeanSpike(DATA);         end     elseif sum(strcmp(type,{'xcorradjacent'}))       AllV.SpikeDraw(a,b,'xcorradj');       elseif sum(strcmp(type,{'xcorrallprobes'}))       AllV.SpikeDraw(a,b,'xcorrall');       elseif strcmp(type,'rateseq')         if strcmp(DATA.plot.expttype,'trialcounts')             DATA.plot.expttype = 'means';             DATA.plot.expt = 0;             set(a,'checked', 'off');         else             DATA.plot.expttype = 'trialcounts'             DATA.plot.expt = 1;             set(a,'checked', 'on');             AllV.PlotExptCounts(DATA);         end            end     set(DATA.toplevel,'UserData',DATA);     