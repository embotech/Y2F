function param = newAdditiveQcqpParam( maps2index, origin, maps2mat, factor )
%NEWQCQPPARAM Creates new additive parameter for QCQPs
% structure:
%     param.maps2index      index of QCQP matrix elements that is affected
%                           by parameter
%     param.origin          cell array filled with two-element vectors. The
%                           first element of each vector contains a
%                           parameter id, the second one the exponent for
%                           this parameter
%     param.factor          factor by which parameter value has to be
%                           multiplied before it is added
%     param.maps2mat        index of QCQP matrix (only relevant for quad.
%                           constraints) 
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech GmbH, Zurich, Switzerland, 2013-2016.

param.maps2index = maps2index;

if nargin >= 2
    param.origin = origin;
else
    param.origin = {};   
end

if nargin >= 3
    param.maps2mat = maps2mat;
else
    param.maps2mat = 1;
end

if nargin >= 4
    param.factor = factor;
else
    param.factor = 1;
end


end

