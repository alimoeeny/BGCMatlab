classdef fullv
    properties (Constant)
        Version = '1.1'
    end
    
    
    
    methods (Static)
            errs = SaveErrors(FullV, varargin);
 [result, FullV] = Check(FullV, varargin)
            errs = ShowErrors(Errors, varargin);
            FullV = ClearErrors(FullV, type, varargin);
                   save(FullV, varargin)
                   split(FullV, varargin)
           FullV = fix(FullV, varargin)

    end
end