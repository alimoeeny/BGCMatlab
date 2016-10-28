function value = objget(obj, val)
%value = objget(obj, val) Gets field val from class properties
f = fields(obj);
if sum(strcmp(val,f))
     value = obj.(val);
else
    value = [];
end