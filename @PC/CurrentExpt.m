function e = CurrentExpt(DATA,varargin)
%CurrentExpt(DATA) return currently selected expt no
%CurrentExpt(DATA,'row') returns row index of expt in cellim

findrow = 0;
 j = 1;
 while j <= length(varargin)
     if strcmp(varargin{j},'row')
         findrow = 1;
     end
     j = j+1;
 end
 
e = DATA.currentpoint(1);
if findrow
    e = find(DATA.CellDetails.exptids == e,1);
end
