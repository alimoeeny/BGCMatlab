function dv = GetdVdY(DATA, varargin)

V = AllV.mygetappdata(DATA,'AllVoltages');
arrstr = GetArrayConfig(DATA.ArrayConfig);
if strcmp(arrstr,'12x2')
    newo = [1:2:24 2:2:24];
    x = diff(V(newo,:,:),1,1);
    dv = x;
    dv(1:2:24,:,:) = x(1:12,:,:);
    dv(2:2:22,:,:) = x(13:end,:,:);
else
    dv = diff(V,1,1);
end

