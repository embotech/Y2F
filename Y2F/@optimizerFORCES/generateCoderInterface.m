function [ success ] = generateCoderInterface( self )
%GENERATEMEXINTERFACECODE generates MEX C code that will prepare the
%parameters for the FORCESPRO solver. It also assembles the correct outputs.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Embotech AG, Zurich, Switzerland, 2013-2023.
    success = 1;
    disp('Generating Coder MATLAB function...');
    success = generateCoderMatlabFunction(self) & success;
    disp('Generating Coder Buildable class...');
    success = generateCoderBuildableClass(self) & success;
    disp('Generating Coder Simulink Block creator script...');
    success = generateCoderSimulinkBlockCreator(self) & success;
    if success
        solverName = self.default_codeoptions.name;
        copyfile([solverName,filesep,'interface',filesep,solverName,'Buildable.m'], '.');
    end
end

function [inputs, outputs] = checkCoderInputsOutputs(self)
    inputs = cell(1, self.numParams);
    for i = 1:self.numParams
        inputs{i} = struct();
        inputs{i}.name = self.paramNames{i};
        if self.paramSizes(i,1) == 1 && self.paramSizes(i,2) == 1 % scalar parameter
            inputs{i}.help_message = 'scalar';
        elseif self.paramSizes(i,1) == 1 % row vector
            inputs{i}.help_message = sprintf('row vector of length %u', self.paramSizes(i,2));
        elseif self.paramSizes(i,2) == 1 % column vector
            inputs{i}.help_message = sprintf('column vector of length %u',self.paramSizes(i,1));
        else
            inputs{i}.help_message = sprintf('matrix of size [%ux%u]',self.paramSizes(i,1),self.paramSizes(i,2));
        end
    end
    inputs = cell2mat(inputs);

    numExtraOutputs = 0;
    if isfield(self.default_codeoptions, 'showinfo') && self.default_codeoptions.showinfo > 0
        numExtraOutputs = 4;
    end

    numOutputs = numel(self.outputSize);
    outputs = cell(1, numOutputs + numExtraOutputs);
    for i = 1:numOutputs
        outputs{i} = struct();
        outputs{i}.name = self.outputNames{i};
        if self.outputSize{i}(1) == 1 && self.outputSize{i}(2) == 1 % scalar output
            outputs{i}.help_message = 'scalar';
        elseif self.outputSize{i}(1) == 1 % row vector
            outputs{i}.help_message = sprintf('row vector of length %u', self.outputSize{i}(2));
        elseif self.outputSize{i}(2) == 1 % column vector
            outputs{i}.help_message = sprintf('column vector of length %u',self.outputSize{i}(1));
        else
            outputs{i}.help_message = sprintf('matrix of size [%ux%u]',self.outputSize{i}(1),self.outputSize{i}(2));
        end
        outputs{i}.copy = sprintf('%s = output.%s;', self.outputNames{i}, self.outputNames{i});
    end
    if isfield(self.default_codeoptions, 'showinfo') && self.default_codeoptions.showinfo > 0
        outputs{numOutputs + 1} = struct('name', 'exitflag', 'help_message', 'scalar: (1) Optimal solution has been found (subject to desired accuracy); (0) Maximum number of interior point iterations reached; (-7) Line search could not progress', 'copy', 'exitflag = min(double(solver_exitflag));');
        outputs{numOutputs + 2} = struct('name', 'iterations', 'help_message', 'iteration number', 'copy', 'iterations = sum(double(info.it));');
        outputs{numOutputs + 3} = struct('name', 'solve_time', 'help_message', 'total solve time', 'copy', 'solve_time = sum(double(info.solvetime));');
        outputs{numOutputs + 4} = struct('name', 'primal_obj', 'help_message', 'primal objective', 'copy', 'primal_obj = sum(double(info.pobj));');
    end
    outputs = cell2mat(outputs);
end

function [inputs, outputs] = getSolverInputsOutputs(self)
    inputs = cell(1, self.numParams);
    for i = 1:self.numParams
        inputs{i} = struct();
        inputs{i}.name = self.paramNames{i};
        inputs{i}.rows = self.paramSizes(i,1);
        inputs{i}.cols = self.paramSizes(i,2);
    end

    outputs = cell(1, numel(self.outputSize));
    for i = 1:numel(self.outputSize)
        outputs{i} = struct();
        outputs{i}.name = self.outputNames{i};
        outputs{i}.rows = self.outputSize{i}(1);
        outputs{i}.cols = self.outputSize{i}(2);
    end
end

function info_fields = getInfoFields(self)
    info_fields = cell(0);
    info_fields{end+1} = struct('name', 'it', 'cast', 'int32');
    info_fields{end+1} = struct('name', 'it2opt', 'cast', 'int32');
    info_fields{end+1} = struct('name', 'res_eq', 'cast', 'double');
    if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
        info_fields{end+1} = struct('name', 'res_ineq', 'cast', 'double');
    else % ADMM has slightly different fields
        info_fields{end+1} = struct('name', 'res_dual', 'cast', 'double');
    end
    info_fields{end+1} = struct('name', 'pobj', 'cast', 'double');
    info_fields{end+1} = struct('name', 'dobj', 'cast', 'double');
    info_fields{end+1} = struct('name', 'dgap', 'cast', 'double');
    info_fields{end+1} = struct('name', 'rdgap', 'cast', 'double');
    if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
        info_fields{end+1} = struct('name', 'mu', 'cast', 'double');
        info_fields{end+1} = struct('name', 'mu_aff', 'cast', 'double');
        info_fields{end+1} = struct('name', 'sigma', 'cast', 'double');
        info_fields{end+1} = struct('name', 'lsit_aff', 'cast', 'int32');
        info_fields{end+1} = struct('name', 'lsit_cc', 'cast', 'int32');
        info_fields{end+1} = struct('name', 'step_aff', 'cast', 'double');
        info_fields{end+1} = struct('name', 'step_cc', 'cast', 'double');
    end
    info_fields{end+1} = struct('name', 'solvetime', 'cast', 'double');
end

function [ success ] = generateCoderMatlabFunction(self)

    % Get solver name from option
    solverName = self.default_codeoptions.name;
    [inputs, outputs] = checkCoderInputsOutputs(self);
    
    exitflag_var = 'solver_exitflag';
    info_struct = 'info';
    if ~isfield(self.default_codeoptions, 'showinfo') || self.default_codeoptions.showinfo <= 0
        exitflag_var = '~';
        info_struct = '~';
    end
    
    fileID = fopen([solverName,filesep,'interface',filesep,solverName,'_coderFunction.m'],'w');

    fprintf(fileID, '%% [OUTPUTS] = %s(INPUTS) solves an optimization problem where:\n', upper(solverName));
    fprintf(fileID, '%% Inputs:\n');
    for i = 1:length(inputs)
        fprintf(fileID, '%% - %s - %s\n', inputs(i).name, inputs(i).help_message);
    end
    fprintf(fileID, '%% Outputs:\n');
    for i = 1:length(outputs)
        fprintf(fileID, '%% - %s - %s\n', outputs(i).name, outputs(i).help_message);
    end
    fprintf(fileID, 'function [%s] = %s(%s)\n', strjoin({outputs.name}, ', '), solverName, strjoin({inputs.name}, ', '));
    fprintf(fileID, '    [output, %s, %s] = %sBuildable.forcesCall(%s);\n', exitflag_var, info_struct, solverName, strjoin({inputs.name}, ', '));
    for i = 1:length(outputs)
        fprintf(fileID, '    %s\n', outputs(i).copy);
    end
    fprintf(fileID, 'end\n');

    fclose(fileID);

    success = 1;
end

function [ success ] = generateCoderBuildableClass(self)

    % Get solver name from option
    solverName = self.default_codeoptions.name;
    [inputs, outputs] = getSolverInputsOutputs(self);
    info_fields = getInfoFields(self);
    
    main_separator = ['...\n', repmat(' ', 1, 12)];
    info_separator = [',...\n', repmat(' ', 1, 26)];
    input_output_separator = [',...\n', repmat(' ', 1, 28)];
    
    exitflag_c_param = 'coder.wref(exitflag_c)';
    if self.numSolvers > 1
        exitflag_c_param = 'exitflag_c';
    end
    fileID = fopen([solverName,filesep,'interface',filesep,solverName,'Buildable.m'],'w');

    fprintf(fileID, 'classdef %sBuildable < coder.ExternalDependency\n', solverName);
    fprintf(fileID, '\n');
    fprintf(fileID, '    methods (Static)\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function name = getDescriptiveName(~)\n');
    fprintf(fileID, '            name = mfilename;\n');
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function b = isSupportedContext(context)\n');
    fprintf(fileID, '            b = context.isMatlabHostTarget();\n');
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function updateBuildInfo(buildInfo, cfg)\n');
    fprintf(fileID, '            buildablepath = fileparts(mfilename(''fullpath''));\n');
    fprintf(fileID, '            [solverpath, foldername] = fileparts(buildablepath);\n');
    fprintf(fileID, '            [~, solvername] = fileparts(solverpath);\n');
    fprintf(fileID, '            %% if the folder structure does not match to the interface folder, we assume it''s the directory that contains the solver\n');
    fprintf(fileID, '            if(~strcmp(foldername, ''interface'') || ~strcmp(solvername, ''%s''))\n', solverName);
    fprintf(fileID, '                solverpath = fullfile(buildablepath, ''%s'');\n', solverName);
    fprintf(fileID, '            end\n');
    fprintf(fileID, '            solverInfo = struct();\n');
    fprintf(fileID, '            solverInfo.solvername = ''%s'';\n', solverName);
    fprintf(fileID, '            solverInfo.solverpath = solverpath;\n');
    fprintf(fileID, '            solverInfo.isY2F = true;\n');
    fprintf(fileID, '            solverInfo.numY2Fsolvers = %d;\n', self.numSolvers);
    fprintf(fileID, '            ForcesUpdateBuildInfo(buildInfo, cfg, solverInfo);\n');
    fprintf(fileID, '            postUpdateBuildInfoScript = [solverInfo.solvername, ''PostUpdateBuildInfo''];\n');
    fprintf(fileID, '            if exist(fullfile(buildablepath, [postUpdateBuildInfoScript, ''.m'']), ''file'')\n');
    fprintf(fileID, '                postUpdateBuildInfo = str2func(postUpdateBuildInfoScript);\n');
    fprintf(fileID, '                postUpdateBuildInfo(buildInfo, cfg, solverInfo);\n');
    fprintf(fileID, '            end\n');
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function [output,exitflag,info] = forcesInitOutputsMatlab()\n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('infos_%s = coder.nullcopy(zeros(1, %d));', s.name, self.numSolvers), info_fields, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            info = struct(%s);\n', strjoin(cellfun(@(s) sprintf('''%s'', infos_%s', s.name, s.name), info_fields, 'UniformOutput', false), info_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('outputs_%s = coder.nullcopy(zeros(%d, %d));', s.name, s.rows, s.cols), outputs, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            output = struct(%s);\n', strjoin(cellfun(@(s) sprintf('''%s'', outputs_%s', s.name, s.name), outputs, 'UniformOutput', false), input_output_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            exitflag = coder.nullcopy(zeros(1,%d));\n', self.numSolvers);
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function [output,exitflag,info] = forcesCallWithParams(params)\n');
    fprintf(fileID, '            [output,exitflag,info] = %sBuildable.forcesCall(%s);\n', solverName, strjoin(cellfun(@(s) sprintf('params.%s', s.name), inputs, 'UniformOutput', false), ', '));
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        \n');
    fprintf(fileID, '        function [output,exitflag,info] = forcesCall(%s)\n', strjoin(cellfun(@(s) s.name, inputs, 'UniformOutput', false), ', '));
    fprintf(fileID, '            solvername = ''%s'';\n', solverName);
    fprintf(fileID, '            \n');
    fprintf(fileID, '            params = struct(%s);\n', strjoin(cellfun(@(s) sprintf('''%s'', double(%s)', s.name, s.name), inputs, 'UniformOutput', false), input_output_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            [output_c, exitflag_c, info_c] = %sBuildable.forcesInitOutputsC();\n', solverName);
    fprintf(fileID, '            \n');
    fprintf(fileID, '            headerName = [solvername ''.h''];\n');
    fprintf(fileID, '            coder.cinclude(headerName);\n');
    fprintf(fileID, '            %% define solver input information (params, file and CasADi)\n');
    fprintf(fileID, '            coder.cstructname(params, [solvername ''_params''], ''extern'', ''HeaderFile'', headerName);\n');
    fprintf(fileID, '            fp = coder.opaque(''FILE *'', ''NULL'', ''HeaderFile'', headerName);\n');
    fprintf(fileID, '            %% define solver output information (output, exitflag, info)\n');
    fprintf(fileID, '            coder.cstructname(output_c,[solvername ''_output''], ''extern'', ''HeaderFile'', headerName);\n');
    fprintf(fileID, '            coder.cstructname(info_c,[solvername ''_info''], ''extern'', ''HeaderFile'', headerName);\n');
    fprintf(fileID, '            coder.ceval([solvername ''_solve''], coder.rref(params), coder.wref(output_c), ...\n');
    fprintf(fileID, '                        %s, coder.wref(info_c), fp);\n', exitflag_c_param);
    fprintf(fileID, '            \n');
    fprintf(fileID, '            [output, exitflag, info] = %sBuildable.forcesInitOutputsMatlab();\n', solverName);
    fprintf(fileID, '            \n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('info.%s = cast(info_c.%s, ''like'', info.%s);', s.name, s.name, s.name), info_fields, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('output.%s = cast(output_c.%s, ''like'', output.%s);', s.name, s.name, s.name), outputs, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            exitflag = cast(exitflag_c, ''like'', exitflag);\n');
    fprintf(fileID, '        end\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    \n');
    fprintf(fileID, '    methods (Static, Access = private)\n');
    fprintf(fileID, '        function [output,exitflag,info] = forcesInitOutputsC()\n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('infos_%s = coder.nullcopy(%s(zeros(1, %d)));', s.name, s.cast, self.numSolvers), info_fields, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            info = struct(%s);\n', strjoin(cellfun(@(s) sprintf('''%s'', infos_%s', s.name, s.name), info_fields, 'UniformOutput', false), info_separator));
    fprintf(fileID, '            \n');
    fprintf(fileID, '            %s\n', strjoin(cellfun(@(s) sprintf('outputs_%s = coder.nullcopy(double(zeros(%d, %d)));', s.name, s.rows, s.cols), outputs, 'UniformOutput', false), main_separator));
    fprintf(fileID, '            output = struct(%s);\n', strjoin(cellfun(@(s) sprintf('''%s'', outputs_%s', s.name, s.name), outputs, 'UniformOutput', false), input_output_separator));
    fprintf(fileID, '            exitflag = coder.nullcopy(int32(zeros(1,%d)));\n', self.numSolvers);
    fprintf(fileID, '        end\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, 'end\n');

    fclose(fileID);

    success = 1;
end

function [ success ] = generateCoderSimulinkBlockCreator(self)

    % Get solver name from option
    solverName = self.default_codeoptions.name;

    [inputs, outputs] = checkCoderInputsOutputs(self);
    
    min_blockheight = 180;
    per_port_blockheight = 17;
    blockheight = min_blockheight + (max(length(inputs), length(outputs)) * per_port_blockheight);
    param_separator = [',...\n', repmat(' ', 1, 29)];
    param_sizes = cell(1, self.numParams);
    for i = 1:self.numParams
        param_sizes{i} = sprintf('''%s'', ''[%d,%d]''', self.paramNames{i}, self.paramSizes(i,1), self.paramSizes(i,2));
    end
    
    fileID = fopen([solverName,filesep,'interface',filesep,solverName,'_createCoderMatlabFunction.m'],'w');

    fprintf(fileID, '%% %s_CREATECODERMATLABFUNCTION(MODELNAME, BLOCKNAME)\n', upper(solverName));
    fprintf(fileID, '%%\n');
    fprintf(fileID, '%% This function generates a Simulink Model with Fixed-step solver type named MODELNAME to add the %s FORCESPRO solver as \n', solverName);
    fprintf(fileID, '%% a MATLAB Function Simulink Block named BLOCKNAME so that it can be used for Simulink simulation and \n');
    fprintf(fileID, '%% Simulink Code Generation.\n');
    fprintf(fileID, '%%\n');
    fprintf(fileID, '%% %s_CREATECODERMATLABFUNCTION(MODELNAME, BLOCKNAME, USEEXISTINGMODEL)\n', upper(solverName));
    fprintf(fileID, '%% can be used if the selected Simulink Model already exists in order to add the Simulink Block to it.\n');
    fprintf(fileID, '%%\n');
    fprintf(fileID, '%% %s_CREATECODERMATLABFUNCTION(MODELNAME, BLOCKNAME, USEEXISTINGMODEL, OPENSIMULINKMODEL)\n', upper(solverName));
    fprintf(fileID, '%% can be used to select whether to open the simulink model.\n');
    fprintf(fileID, 'function %s_createCoderMatlabFunction(modelname, blockname, useExistingModel, openSimulinkModel)\n', solverName);
    fprintf(fileID, '\n');
    fprintf(fileID, '    if nargin < 1 || isempty(modelname)\n');
    fprintf(fileID, '        modelname = ''%s_model'';\n', solverName);
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if ~isvarname(modelname)\n');
    fprintf(fileID, '        error(''Modelname must be a valid variable name'');\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if nargin < 2 || isempty(blockname)\n');
    fprintf(fileID, '        blockname = ''%s_block'';\n', solverName);
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if ~isvarname(blockname)\n');
    fprintf(fileID, '        error(''Blockname must be a valid variable name'');\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if nargin < 3\n');
    fprintf(fileID, '        useExistingModel = false;\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if ~islogical(useExistingModel)\n');
    fprintf(fileID, '        error(''useExistingModel must be a bool'');\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if nargin < 4\n');
    fprintf(fileID, '        openSimulinkModel = true;\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if ~islogical(openSimulinkModel)\n');
    fprintf(fileID, '        error(''openSimulinkModel must be a bool'');\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    \n');
    fprintf(fileID, '    function_name = ''%s_coderFunction.m'';\n', solverName);
    fprintf(fileID, '    position = [170, 99, 650, %d];\n', blockheight);
    fprintf(fileID, '    parameter_sizes = struct(%s);\n', strjoin(param_sizes, param_separator));
    fprintf(fileID, '    \n');
    fprintf(fileID, '    result = exist(modelname, ''file'');\n');
    fprintf(fileID, '    if result ~= 4 && result ~= 0\n');
    fprintf(fileID, '        error(''%%s exists but is not a Simulink model. Please use a different name.'', modelname);\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if result == 0\n');
    fprintf(fileID, '        new_system(modelname);\n');
    fprintf(fileID, '        set_param(modelname, ''SolverType'', ''Fixed-step'')\n');
    fprintf(fileID, '    elseif ~useExistingModel\n');
    fprintf(fileID, '        error(''Simulink Model %%s already exists. Please call the script again with the parameter useExistingModel set to true to add MATLAB Function Simulink Block to existing model.'', modelname);\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if openSimulinkModel\n');
    fprintf(fileID, '        open_system(modelname);\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '\n');
    fprintf(fileID, '    blockpath = [modelname, ''/'', blockname];\n');
    fprintf(fileID, '    blockExists = true;\n');
    fprintf(fileID, '    try\n');
    fprintf(fileID, '        get_param(blockpath, ''ObjectParameters'');\n');
    fprintf(fileID, '    catch\n');
    fprintf(fileID, '        blockExists = false;\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    if blockExists\n');
    fprintf(fileID, '        error(''Blockname %%s already exists in Simulink Model %%s. Please select a different blockname'', blockname, modelname);\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    add_block(''simulink/User-Defined Functions/MATLAB Function'', blockpath);\n');
    fprintf(fileID, '    set_param(blockpath, ''Position'', position);\n');
    fprintf(fileID, '\n');
    fprintf(fileID, '    cur_dir = fileparts(mfilename(''fullpath''));\n');
    fprintf(fileID, '    script_filepath = fullfile(cur_dir, function_name);\n');
    fprintf(fileID, '    try\n');
    fprintf(fileID, '        blockconfig = find(Simulink.Root,''-isa'',''Stateflow.EMChart'', ''Path'', blockpath);\n');
    fprintf(fileID, '        blockconfig.Script = fileread(script_filepath);\n');
    fprintf(fileID, '        parameters = fieldnames(parameter_sizes);\n');
    fprintf(fileID, '        for i = 1:length(parameters)\n');
    fprintf(fileID, '            inputportconfig = find(Simulink.Root,''-isa'',''Stateflow.Data'', ''Path'', blockpath, ''Name'', parameters{i});\n');
    fprintf(fileID, '            inputportconfig.Props.Array.Size = parameter_sizes.(parameters{i});\n');
    fprintf(fileID, '        end\n');
    fprintf(fileID, '        blockconfig.Locked = 1;\n');
    fprintf(fileID, '    catch\n');
    fprintf(fileID, '        error(''This Simulink version does not support programmatic editing of the script of the MATLAB Function Simulink Block. Please manually create a MATLAB Function block, copy the script from ''''%%s'''' and set the following dimensions for the input ports:\\n%%s.'', script_filepath, printCharStruct(parameter_sizes));\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, '    save_system(modelname);\n');
    fprintf(fileID, 'end\n');
    fprintf(fileID, '\n');
    fprintf(fileID, 'function message = printCharStruct(struct_var)\n');
    fprintf(fileID, '    message = '''';\n');
    fprintf(fileID, '    fields = fieldnames(struct_var);\n');
    fprintf(fileID, '    for i = 1:length(fields)\n');
    fprintf(fileID, '        message = sprintf(''%%s''''%%s'''': ''''%%s''''\\n'', message, fields{i}, struct_var.(fields{i}));\n');
    fprintf(fileID, '    end\n');
    fprintf(fileID, 'end\n');

    fclose(fileID);

    success = 1;
end
