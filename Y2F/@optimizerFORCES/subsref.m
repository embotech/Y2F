function [output,exitflag,info] = subsref(self,subs)
%SUBSREF Overload of subsref. Makes A{B} for object A possible.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech GmbH, Zurich, Switzerland, 2013-2016.

if isequal(subs.type,'()')
    error('Cannot index in OPTIMIZERFORCES. Perhaps you mean {}');
elseif isequal(subs.type,'.')
    error('No fields accessible in OPTIMIZERFORCES.');
elseif isequal(subs.type,'{}') % --> call solver
    if ~isa(subs.subs{1},'cell')
        paramValues = subs.subs;
    else
        paramValues = subs.subs{1};
    end
    [output,exitflag,info] = self.interfaceFunction(paramValues);
end