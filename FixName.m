function name = FixName(name, varargin)

os = computer;
if strncmp(os,'PCW',2)
else
   name = strrep(name,'\','/') ;
end