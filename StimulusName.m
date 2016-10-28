function name = StimulusName(st)
%Convert stimulus name string to num or vice versa
%binoc.StimulusName is the new version of this
stimnames = {'none' ,	'gabor',	'rds' ,	'grating',	'bar',	'circle',...
	'rectangle','test',	'square',	  'probe',	  '2grating',  'cylinder',...
	  'corrug',	'sqcorrug',	'twobar',	'rls', 'annulus', 'rdssine', 'nsines', 'rlssine',...
	  'radial', 'image', 'checker'};
  
if ischar(st)
    name = find(strcmp(st, stimnames));
    if ~isempty(name)
        name = name-1;
    else
        name = NaN;
    end
elseif isfield(st,'Header') % and expt strcut
    name = StimulusName(GetEval(st,'st'));
    if strcmp(name,'unknown') && length(st.Trials) > 10
        for j = 1:length(stimnames)
            if ~isempty(strfind(st.Header.name,stimnames{j}))
                name = stimnames{j};
            end
        end
    end
elseif isnan(st)
    name = 'unknown';
else
    name = stimnames{st+1};
end
