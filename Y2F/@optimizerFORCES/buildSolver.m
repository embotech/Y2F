function [ success ] = buildSolver( self )
%BUILDSOLVER Generates FORCESPRO solver(s) as well as MEX wrapper to call
%main controller.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech AG, Zurich, Switzerland, 2013-2023.

    % ensure minimum FORCESPRO version
    [version, ~] = FORCESversion();
    reqVersion = '6.2.0';
    if (~satisfiesMinimalVersion(version, reqVersion))
        error(['Y2F requires FORCESPRO v', reqVersion ,' or higher (installed version: v', version, ')']);
    end

    success = 1;
    for i=1:self.numSolvers
        success = generateCode( self.stages{i},self.params{i},self.codeoptions{i},self.outputFORCES{i} ) & success;
    end
    if ~success
        error('Code generation was not successful');
    end

    %% Generate MEX code that is called when the solver is used
    success = generateY2FInterfaces( self ) & success;

end


function [ success ] = generateY2FInterfaces( self )
% Helper function to generate MEX code that is called when the solver is used.

    success = 1;

    disp('Generating C interface...');
    %generateSolverInterfaceCode(sys);
    success = generateCInterfaceCode(self) & success;
    success = generateMEXInterfaceCode(self) & success;
    success = generateSimulinkInterfaceCode(self) & success;
    success = packageSolverCode(self) & success;

    % Generate help file
    disp('Writing help file...');
    success = generateHelp(self) & success;

    % Compile MEX code
    disp('Compiling MEX code for solver interface...');
    success = compileMEXInterfaceCode(self, '_mex', '') & success;

    if (~isfield(self.default_codeoptions,'BuildSimulinkBlock') || self.default_codeoptions.BuildSimulinkBlock ~= 0)
        % Compile Simulink code (is optional)
        disp('Compiling Simulink code for solver interface...');
        success = compileMEXInterfaceCode(self, '_simulinkBlock', '_simulinkBlock') & success;

        % Compile Simulink code
        disp('Generating Simulink Block...');
        success = generateSimulinkBlock(self) & success;
    end
    
    % Generate coder interface files
    disp('Generating coder interface files...');
    success = generateCoderInterface(self) & success;

end
