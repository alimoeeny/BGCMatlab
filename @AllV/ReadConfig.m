function DATA = ReadConfig(DATA, name, varargin)
%AllV.ReadConfig. Read Config and then do a couple of AllV Specific checks

auto = DATA.auto;
DATA = ReadConfig(DATA, name, varargin{:}); %%NOT recursive. Calls verion in ~/matlab
DATA.auto.loadfromspikes = auto.loadfromspikes;
if DATA.interactive < 0
    DATA.auto.uselastcluster = 0;
end