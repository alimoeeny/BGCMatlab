function DATA = SetSpaceFromAxis(DATA, ax, varargin)
%AllV.SetSpaceFromAxis use axis user data to define cluster space

pcplot = get(gca,'UserData');
if isfield(pcplot,'Variables')
    DATA.elmousept.pcplot = pcplot.pcplot;
    DATA.elmousept.Variables = pcplot.Variables;
else
    DATA.elmousept.pcplot = pcplot;
    DATA.elmousept.Variables = {};
end