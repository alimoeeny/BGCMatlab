 function SaveExtras(DATA)     Xbits = [];     outname = [DATA.name '/PlotClusterExtra.mat'];     if isfield(DATA,'GaussFitdp')          Xbits.GaussFitdp = DATA.GaussFitdp;         Xbits.gmfitpos = DATA.gmfitpos;     end     if isfield(DATA,'mahal3')         Xbits.mahal3 = DATA.mahal3;     end     if isfield(DATA,'xcorrs')         Xbits.xcorrs= DATA.xcorrs;         Xbits.xcorrval = DATA.xcorrval;     end     if isfield(DATA,'xysdindex')         Xbits.xysdindex = DATA.xysdindex;     end     names = {'shiftmatrix' 'fitjumps'};     for j = 1:length(names)     if isfield(DATA,names{j})         Xbits.(names{j}) = DATA.(names{j});     end     end     if ~isempty(Xbits)         save(outname,'-v7.3','Xbits')     end