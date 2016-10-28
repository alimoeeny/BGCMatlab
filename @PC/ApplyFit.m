function [C ] = ApplyFit(DATA, C, fitnumber, varargin )

res = PC.CallAllVPcs(DATA, C.exptid,.C.probe(1), 'setfit',fitnumber);


