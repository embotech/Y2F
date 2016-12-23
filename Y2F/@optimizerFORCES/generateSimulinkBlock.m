function generateSimulinkBlock( self )
%GENERATESIMULINKBLOCK Create a model containing a Simulink block for the
%solver
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech GmbH, Zurich, Switzerland, 2013-2016.

success = 0;

% Get solver name from option
solverName = self.default_codeoptions.name;

% Check if FORCES solver has been generated
if ~isdir(solverName)
    error('Solver ''%s'' has not been generated!', solverName)
end

model = [solverName '_lib'];
block = solverName;

% Create Simulink library
if exist(model)
    close_system(model, 0);
end
new_system(model, 'Library');

% Add S-function block
add_block('built-in/S-Function', [model '/' block]);
set_param([model '/' block], 'FunctionName', [solverName '_simulinkBlock']);

% Create mask to name input/output ports and add image
mask = Simulink.Mask.create([model '/' block]);

iconDrawingString = 'image(''forcesprologo.jpg'', ''center'', ''on'')';
for i=1:self.numParams
    iconDrawingString = sprintf('%s;port_label(''input'', %u, ''%s'')', iconDrawingString, i, self.paramNames{i});
end
for i=1:numel(self.outputSize)
    iconDrawingString = sprintf('%s;port_label(''output'', %u, ''%s'')', iconDrawingString, i, self.outputNames{i});
end
set_param([model '/' block], 'MaskDisplay', iconDrawingString)

% Set position of block
set_param([model '/' block], 'Position', [170, 99, 550, 200])

% Generate description and help
desc = sprintf(['---- Simulink block encapsulating your customized solver %s ----\n\n' ...
    '%s : A fast customized optimization solver.'], solverName, solverName);
set_param([model '/' block], 'MaskDescription', desc);

help = sprintf('%s_simulinkBlock provides an easy Simulink interface for simulating your customized solver.\n\n',solverName);
help = sprintf('%sOUTPUTS = %s(INPUTS) solves an optimization problem where:\n\n', help, solverName);
help = sprintf('%sINPUTS:\n', help);

for i=1:self.numParams
    if self.paramSizes(i,1) == 1 && self.paramSizes(i,2) == 1 % scalar parameter
        help = sprintf('%s- %s (a scalar)\n',help,self.paramNames{i});
    elseif self.paramSizes(i,1) == 1 % row vector
        help = sprintf('%s- %s (row vector of length %u)\n',help,self.paramNames{i},self.paramSizes(i,2));
    elseif self.paramSizes(i,2) == 1 % column vector
        help = sprintf('%s- %s (column vector of length %u)\n',help,self.paramNames{i},self.paramSizes(i,1));
    else
        help = sprintf('%s- %s (matrix of size [%u x %u])\n',help,self.paramNames{i},self.paramSizes(i,1),self.paramSizes(i,2));
    end
end

help = sprintf('%s\nOUTPUTS:\n',help);

for i=1:numel(self.outputSize)
    if self.outputSize{i}(1) == 1 && self.outputSize{i}(2) == 1 % scalar parameter
        help = sprintf('%s- %s (a scalar)\n',help,self.outputNames{i});
    elseif self.outputSize{i}(1) == 1 % row vector
        help = sprintf('%s- %s (row vector of length %u)\n',help,self.outputNames{i},self.outputSize{i}(2));
    elseif self.outputSize{i}(2) == 1 % column vector
        help = sprintf('%s- %s (column vector of length %u)\n',help,self.outputNames{i},self.outputSize{i}(1));
    else
        help = sprintf('%s- %s (matrix of size [%u x %u])\n',help,self.outputNames{i},self.outputSize{i}(1),self.outputSize{i}(2));
    end
end

help = sprintf('%s\nFor more information, see https://www.embotech.com/FORCES-Pro/How-to-use/Simulink-Interface/Simulink-Block',help);

set_param([model '/' block], 'MaskHelp', help);

% Save system
save_system(model);

end

