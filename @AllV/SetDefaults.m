function DATA = SetDefaults(DATA, defmode, varargin)
%DATA = AllV.SetDefaults(DATA, defmode, varargin)
%         DATA, 'safeauto') turns off options like uselastcluster
%              that must be set from the GUI - will mess up reclassify
%              otherwise
%     



if nargin < 2
    defmode = 'safeauto';
end

if strcmp(defmode,'safeauto')
    DATA.fullvswitchmode.applylast = 0;
    DATA.auto.uselastcluster = 0;
    DATA.auto.advanceexpt = 0;
end