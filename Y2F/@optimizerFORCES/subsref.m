function [output,exitflag,info] = subsref(self,subs)
%SUBSREF Overload of subsref. Makes A{B} for object A possible.

if isequal(subs.type,'()')
    error('Cannot index in OPTIMIZERFORCES. Perhaps you mean {}');
elseif isequal(subs.type,'.')
    error('No fields accessible in OPTIMIZERFORCES.');
elseif isequal(subs.type,'{}') % --> call solver
    if ~isa(subs.subs{1},'cell')
        paramValues = {subs.subs{1}};
    else
        paramValues = subs.subs{1};
    end
    [output,exitflag,info] = self.interfaceFunction(paramValues);
end