function [sys, success] = optimizerFORCES( constraint,objective,codeoptions,parameters,solverOutputs,parameterNames,outputNames )
%OPTIMIZERFORCES Generates a FORCES Pro solver from a YALMIP problem formulation
%
%   solver = OPTIMIZERFORCES(constraint,objective,codeoptions,parameters,solverOutputs)
%   generates a solver using FORCES Pro that solves the specified problem
%   for given parameter values and returns the value of the specified
%   outputs. The solver can be called  by indexing the returned OPTIMIZERFORCES
%   object with the parameter values. The parameter values can be passed in
%   a number of different ways:
%   
%       output = solver{param1, param2, ...} OR
%       output = solver{ {param1, param2, ...} } OR
%       output = solver(param1, param2, ...)
%   
%   This API is compatible with the YALMIP function OPTIMIZER. To obtain more
%   information on the solver, you can type 'help <solvername>' after the code
%   has been generated, where <solvername> is the name you give to the solver
%   via codeoptions.name.
%
%
%   solver = OPTIMIZERFORCES(constraint,objective,codeoptions,parameters,solverOutputs,parameterNames,outputNames)
%   same as above, but the generated code and all its associated help
%   files, the Simulink block etc. will carry the names provided in the
%   last two arguments for parameters and outputs, respectively. It is
%   strongly recommended to provide these names.
%
%
%   Example usage:
%
%       sdpvar x y p
%       options = getOptions('solver_name')
%     	solver = optimizerFORCES([x + y == p], x^2+y^2, options, p, {x,y}, {'p'}, {'x_optimal','y_optimal'})
%       solution = solver{5} % solve problem for p = 5
%
%
%   Inputs:
%       constraint:     constraints in YALMIP format, see e.g. OPTIMIZER
%       objective:      objective in YALMIP format, see e.g. OPTIMIZER
%       codeoptions:    options for FORCES Pro solver, see GETOPTIONS
%       parameters:     Single SDPVAR object or cell array of SDPVAR
%                       objects that should be considered a parameter
%       solverOutputs:  Single SDPVAR object or cell array of SDPVAR
%                       objects whose value(s) should be returned by the
%                       solver, can be a linear combination of decision
%                       variables and parameters
%       parameterNames: (optional) Cell array of strings with names for
%                       parameters that will, for example, be used in the
%                       Simulink block. If names are not specified, they
%                       will be auto-generated.
%       outputNames:    (optional) Cell array of strings with names for
%                       outputs that will, for example, be used in the
%                       Simulink block. If names are not specified, they
%                       will be auto-generated.
%
%   Outputs:
%       solver:         reference to OPTIMIZERFORCES object. Use this to
%                       call solver. Example:
%                           solver = optimizerFORCES(...);
%                           x = solver{paramValues};
%       success:        1 if solver generation was successful, 0 otherwise
%
%
%   Getting details on solver:
%
%       help <solvername> % get more details on the generated solver
%
% See also SDPVAR OPTIMIZER
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech GmbH, Zurich, Switzerland, 2013-2016.

disp('YALMIP-to-FORCES code generator')
disp('-------------------------------')

% Check if YALMIP is installed
if ~exist('optimizer','file')
    error('YALMIP could not be found. Please make sure it is installed correctly.')
end

% We need all arguments
switch nargin
    case 0
        error('Constraints not found')
    case 1
        error('Objective not found')
    case 2
        error('Solver options not found')
    case 3
        error('Parameter(s) not found')
    case 4
        error('Output(s) not found')
end

% Make valid solver name
if ~verLessThan('matlab', '8.3')
    codeoptions.name = matlab.lang.makeValidName(codeoptions.name);
else
    codeoptions.name = genvarname(codeoptions.name);
end

% We need parameters
if isempty(parameters)
    error('FORCES Pro does not support problems without parameters.');
end

% Read parameter names if they were passed along
if nargin >= 6
    if ~iscellstr(parameterNames)
        error('parameterNames needs to be a cell array of strings.')
    end
    
    % Fix names (make them valid and unique)
    if ~verLessThan('matlab', '8.3')
        parameterNames = matlab.lang.makeValidName(parameterNames);
        parameterNames = matlab.lang.makeUniqueStrings(parameterNames);
    else
        parameterNames = genvarname(parameterNames);
    end
elseif isa(solverOutputs,'sdpvar')
    % Single parameter supplied, we might get its name!
    name = inputname(4);
    if ~isempty(name)
        parameterNames = {name};
    else
        parameterNames = {};
        warning('Y2F:noParameterNames',['No parameter names specified for solver. We recommend adding names for better code documentation. ' ...
        'For more info type ''help optimizerFORCES''.']);
    end
else
    parameterNames = {};
    warning('Y2F:noParameterNames',['No parameter names specified for solver. We recommend adding names for better code documentation. ' ...
        'For more info type ''help optimizerFORCES''.']);
end

% Read output names if they were passed along
if nargin >= 7
    if ~iscellstr(parameterNames)
        error('outputNames needs to be a cell array of strings.')
    end
    
    % Fix names (make them valid and unique)
    if ~verLessThan('matlab', '8.3')
        outputNames = matlab.lang.makeValidName(outputNames);
        outputNames = matlab.lang.makeUniqueStrings(outputNames);
    else
        outputNames = genvarname(outputNames);
    end
elseif isa(solverOutputs,'sdpvar')
    % Single output supplied, we might get its name!
    name = inputname(5);
    if ~isempty(name)
        outputNames = {name};
    else
        outputNames = {};
        warning('Y2F:noOutputNames',['No output names specified for solver. We recommend adding names for better code documentation. ' ...
        'For more info type ''help optimizerFORCES''.']);
    end
else
    outputNames = {};
    warning('Y2F:noOutputNames',['No output names specified for solver. We recommend adding names for better code documentation. ' ...
        'For more info type ''help optimizerFORCES''.']);
end

% Create missing parameter names
while numel(parameterNames) < numel(parameters)
    parameterNames{end+1} = sprintf('param%u', numel(parameterNames)+1);
end

% We allow single parameters --> wrap them in a cell array
if ~iscell(parameters) 
    parameters = {parameters}; % put single param into cell array
end

% Prepare struct that is going to be converted into the optimizerFORCES
% class
sys = struct;

sys.outputIsCell = 1;
if ~iscell(solverOutputs)
    solverOutputs = {solverOutputs};
    sys.outputIsCell = 0;
end

% Create missing parameter names
while numel(outputNames) < numel(solverOutputs)
    outputNames{end+1} = sprintf('output%u', numel(outputNames)+1);
end


%% Call YALMIP and convert QP into FORCES format
disp('This is Y2F (v0.1.10), the YALMIP interface of FORCES Pro.');
disp('For more information visit https://github.com/embotech/y2f');
fprintf('\nUsing YALMIP to convert problem into QP...')
tic;
[internalmodel,H,f,Aineq,bineq,Aeq,beq,lb,ub] = getQpAndModelFromYALMIP();
yalmiptime=toc;
fprintf('   [OK, %5.1f sec]\n', yalmiptime);

% Check if matrices are numeric
if ~isnumeric(H) || ~isnumeric(f)
    error('Y2F can only handle numeric inputs. There are non-numeric terms in the cost.')
end
if ~isnumeric(Aineq) || ~isnumeric(bineq)
    error('Y2F can only handle numeric inputs. There are non-numeric terms in the inequality contraints.')
end
if ~isnumeric(Aeq) || ~isnumeric(beq)
    error('Y2F can only handle numeric inputs. There are non-numeric terms in the equality contraints.')
end
if ~isnumeric(lb) || ~isnumeric(ub)
    error('Y2F can only handle numeric inputs. There are non-numeric terms in the bounds.')
end

% Check if matrices are real
if ~isreal(H) || ~isreal(f)
    error('Y2F can only handle real inputs. There are complex terms in the cost.')
end
if ~isreal(Aineq) || ~isreal(bineq)
    error('Y2F can only handle real inputs. There are complex terms in the inequality contraints.')
end
if ~isreal(Aeq) || ~isreal(beq)
    error('Y2F can only handle real inputs. There are complex terms in the equality contraints.')
end
if ~isreal(lb) || ~isreal(ub)
    error('Y2F can only handle real inputs. There are complex terms in the bounds.')
end

% Check if matrices are doubles
if ~isa(H,'double') || ~isa(f,'double')
    warning('Y2F:nonDoubleCost', 'Y2F can only handle inputs of type ''double''. The cost will be cast to ''double''.')
    H = double(H);
    f = double(f);
end
if ~isa(Aineq,'double') || ~isa(bineq,'double')
    warning('Y2F:nonDoubleInequality', 'Y2F can only handle inputs of type ''double''. The inequality constraints will be cast to ''double''.')
    Aineq = double(Aineq);
    bineq = double(bineq);
end
if ~isa(Aeq,'double') || ~isa(beq,'double')
    warning('Y2F:nonDoubleEquality', 'Y2F can only handle inputs of type ''double''. The equality constraints will be cast to ''double''.')
    Aeq = double(Aeq);
    beq = double(beq);
end
if ~isa(lb,'double') || ~isa(ub,'double')
    warning('Y2F:nonDoubleBounds', 'Y2F can only handle inputs of type ''double''. The bounds will be cast to ''double''.')
    lb = double(lb);
    ub = double(ub);
end

%% Assemble parameters & convert quadratic variables
% Quadratic inequalities are not recognized by YALMIP
% Information is stored in internalmodel

fprintf('Extract parameters and quadratic inequalities from YALMIP model...')
tic;
[qcqpParams,Q,l,r,solverVars,paramVars,yalmipParamMap] = buildParamsAndQuadIneqs();
extractStagesTime=toc;
fprintf('   [OK, %5.1f sec]\n', extractStagesTime);

%%

fprintf('Assembling stages...')
tic;
% Construct matrices where parametric elements are == 1
% This is necessary to build graph and recognise infeasible problems
H_temp = H;
H_temp([qcqpParams.H.maps2index]) = 1;
Aineq_temp = Aineq;
Aineq_temp([qcqpParams.Aineq.maps2index]) = 1;
Aeq_temp = Aeq;
Aeq_temp([qcqpParams.Aeq.maps2index]) = 1;
Q_temp = Q;
for i=1:numel(qcqpParams.Q)
    Q_temp{qcqpParams.Q(i).maps2mat}(qcqpParams.Q(i).maps2index) = 1;
end
l_temp = l;
for i=1:numel(qcqpParams.l)
    l_temp(qcqpParams.l(i).maps2index,qcqpParams.l(i).maps2mat) = 1;
end

%% Warn the user if the problem is/might be infeasible

checkQcqpForInfeasibility()


%% Generate standard stages
% Construct (potentially multiple) path graphs from Qp

% Compute decision variable indices that are needed in output (only
% relevant for separable problems)
outputIdx = [];
for i=1:numel(solverOutputs)
    outputVars = getvariables(solverOutputs{i});
    outputIdx = [outputIdx find(ismember(solverVars, outputVars))];
end

graphComponents = pathGraphsFromQcqp(H_temp,Aineq_temp,Aeq_temp,Q_temp,l_temp);
[graphComponents, stages, params, standardParamValues,forcesParamMap] = stagesFromPathGraphs(graphComponents,H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,qcqpParams,yalmipParamMap,outputIdx);

if numel(stages) == 1
    fprintf('Found %u stages.\n', numel(stages{1}));
else  % we have multiple solvers
    fprintf('The problem is separable. %u solvers are needed:\n', numel(stages));
    for i=1:numel(stages)
        fprintf('    - Solver %u has %u stages\n', i, numel(stages{i}));
    end
end


%% Assemble the rest of the FORCES parameters
% Fake a parameter for each solver if there are none (we need one for FORCES)
sys.solverHasParams = zeros(1,numel(stages));
for i=1:numel(stages)
    if isempty(params{i})
        params{i}(1) = newParam('p',1,'cost.f');
        standardParamValues{i} = stages{i}(1).cost.f;
        stages{i}(1).cost.f = [];
    else % count params
        sys.solverHasParams(i) = 1;        
    end
end

% Mark solvers that contain binary variables
sys.solverIsBinary = zeros(1,numel(stages));
if ~isempty(qcqpParams.bidx) 
    for i=1:numel(stages)
        if any(ismember(cell2mat(graphComponents{i}.vertices),qcqpParams.bidx))
            sys.solverIsBinary(i) = 1;      
        end
    end
end

% Assemble outputs
outputFORCES = buildOutput();

assembleStagesTime=toc;
fprintf('   [OK, %5.1f sec]\n', assembleStagesTime);

%% Generate solver using FORCES
disp('Generating solver using FORCES...')

% set flag to let FORCES know that request came from Y2F
codeoptions.interface = 'y2f';

success = 1;
default_codeoptions = codeoptions;
codeoptions = cell(1,numel(stages));
for i=1:numel(stages)
    % new name for each solver
    codeoptions{i} = default_codeoptions;
    codeoptions{i}.name = sprintf('internal_%s_%u',default_codeoptions.name,i);
    codeoptions{i}.nohash = 1; % added by AD to avoid problem - exeprimental
    success = generateCode(stages{i},params{i},codeoptions{i},outputFORCES{i}) & success;
end

if ~success
    error('Code generation was not successful');
end

% Store temporary data in object
sys.stages = stages;
sys.numSolvers = numel(stages);
sys.params = params;
sys.paramNames = parameterNames;
sys.outputNames = outputNames;
sys.outputFORCES = outputFORCES;
sys.qcqpParams = qcqpParams;
sys.standardParamValues = standardParamValues;
sys.solverVars = solverVars;
sys.parameters = parameters;
sys.forcesParamMap = forcesParamMap;
sys.codeoptions = codeoptions;
sys.default_codeoptions = default_codeoptions;
sys.interfaceFunction = str2func(default_codeoptions.name);

sys = class(sys,'optimizerFORCES');

% Generate MEX code that is called when the solver is used
disp('Generating C interface...');
%generateSolverInterfaceCode(sys);
generateCInterfaceCode(sys);
generateMEXInterfaceCode(sys);
generateSimulinkInterfaceCode(sys);

% Compile MEX code
disp('Compiling MEX code for solver interface...');

compileSolverInterfaceCode(sys);

% Generate help file
disp('Writing help file...');

generateHelp(sys);

if (~isfield(default_codeoptions,'BuildSimulinkBlock') || default_codeoptions.BuildSimulinkBlock ~= 0)
    % Compile Simulink code (is optional)
    disp('Compiling Simulink code for solver interface...');

    compileSimulinkInterfaceCode(sys);

    % Compile Simulink code
    disp('Generating Simulink Block...');

    generateSimulinkBlock(sys);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [internalmodel,H,f,Aineq,bineq,Aeq,beq,lb,ub] = getQpAndModelFromYALMIP()
    % Helper function that uses YALMIP to create Qp from user's constraints
    % and objective. Parameters and quadratic constraints are ignored for
    % now and handled later on.
        
        % Call YALMIP's export
        % model contains Qp (in quadprog format)
        % internalmodel contains data that we need to recover parameters and
        % quadratic inequalities
        options = sdpsettings('solver','+quadprog','verbose',2);
        [model,~,diagnostics,internalmodel] = export(constraint,objective,options,[],[],1);
        if ~isempty(diagnostics)
            disp(diagnostics)
            error('YALMIP was not able to convert the problem into a QCQp.')
        end

        % Get matrices from quadprog model
        H = full(model.Q);
        f = full(model.c);
        Aineq = full(model.A);
        bineq = full(model.b);
        Aeq = full(model.Aeq);
        beq = full(model.beq);
        lb = full(model.lb);
        ub = full(model.ub);

        % YALMIP doesn't always recognize bounds as such (and makes them
        % inequalities) --> we have to convert them
        bounds_idx = []; % remember to delete these
        for i=1:size(Aineq,1)
            vars = find(Aineq(i,:)); % variables used in inequality
            if length(vars) == 1 && internalmodel.variabletype(vars) == 0 % a single linear variable
                % we can make a bound out of this
                bounds_idx(end+1) = i;
                if Aineq(i,vars) > 0 % upper bound
                    ub(vars) = min(ub(vars),bineq(i)/Aineq(i,vars));
                else % < 0 --> lower bound
                    lb(vars) = max(lb(vars),bineq(i)/Aineq(i,vars));
                end
            end
        end
        Aineq(bounds_idx,:) = [];
        bineq(bounds_idx) = [];

    end

    function [qcqpParams,Q,l,r,solverVars,paramVars,yalmipParamMap] = buildParamsAndQuadIneqs()
    % Helper function that builds QCQp parameter list and recognises
    % quadratic inequalities
    
        % Build empty additive QCQp param struct
        qcqpParams.H = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.f = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.Aineq = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.bineq = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.Aeq = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.beq = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.l = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.Q = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.r = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.lb = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.ub = struct('maps2index',{},'maps2origparam',{},...
            'maps2mat',{},'factor',{});
        qcqpParams.bidx = [];

        % Build empty quadratic constraints
        Q = {};
        l = [];
        r = [];

        % Find all YALMIP variable indices for parameters & build parameter map
        paramVars = []; % YALMIP variables that are parameters
        yalmipParamMap = zeros(2,0); % 1st row: index of matrix with values,
                                         % 2nd row: index of element inside matrix
        if ~isempty(parameters)
            sys.paramSizes = zeros(numel(parameters),2);
            for i=1:numel(parameters)
                if ~isa(parameters{i}, 'sdpvar')
                    error('Parameters must be a SDPVAR or a cell array of SDPVARs.');
                end
                
                % store size for code generation
                sys.paramSizes(i,:) = size(parameters{i});
                
                % find YALMIP variables that make up parameter
                newParams = getvariables(parameters{i});
                paramVars = [paramVars newParams];

                % find element inside matrix (given by user) that contains value of
                % parameter
                for p=newParams
                    yalmipParamMap(:,end+1) = [i; find(getbasematrix(parameters{i},p),1)];
                end
            end

            assert(length(paramVars) == size(yalmipParamMap,2));
            sys.numParams = numel(parameters);
        else
            sys.numParams = 0;
        end



        quadIneq = zeros(2,0); % data structure for keeping track of (parametric) quadratic inequalities
                               % first row:     id of linear inequality
                               % second row:    id of quad inequality 

        % Go through all variables and classify them. Remove parameters and
        % interpret 
        removeIdx = []; % don't forget to delete these
        for i=length(internalmodel.used_variables):-1:1 % go through all YALMIP variables
            var = recover(internalmodel.used_variables(i)); % get corresponding SDPVAR object

            % is this variables binary?
            binary = ismember(i,internalmodel.integer_variables);
            if binary
                qcqpParams.bidx(end+1) = i;

                if lb(i) ~= 0 || ub(i) ~= 1
                    error('No integer variables other than binary supported.');
                end
            end

            % Is this variable a linear parameter?
            if is(var, 'linear') && ...
                    any(paramVars == internalmodel.used_variables(i)) 

                % Parameters cannot be binary
                if binary
                    error('Parameters cannot be binary.');
                end

                % index in parameter list (needed to recover value later on)
                p_idx = find(paramVars == internalmodel.used_variables(i),1); 

                % Variable is parameter --> remove it
                removeIdx = [removeIdx i];

                % Check if parameter appears in H --> put into f
                rows = find(H(:,i));
                rows = rows(rows ~= i); % param*param just adds a constant term to cost
                for row=rows'
                    if ~any(paramVars == internalmodel.used_variables(row)) % param1*param2 just adds a constant term to cost
                        qcqpParams.f(end+1) = newAdditiveQcqpParam(row,p_idx,1,0.5*H(row,i));
                    end
                end
                cols = find(H(i,:));
                cols = cols(cols ~= i); % param*param just adds a constant term to cost
                for col=cols
                    if ~any(paramVars == internalmodel.used_variables(col)) % param1*param2 just adds a constant term to cost
                        qcqpParams.f(end+1) = newAdditiveQcqpParam(col,p_idx,1,0.5*H(i,col));
                    end
                end

                % We can ignore it if param appears in f (just a constant term)

                % Check if parameter is used in Aineq --> add to bineq
                rows = find(Aineq(:,i));
                for row=rows'
                    qcqpParams.bineq(end+1) = newAdditiveQcqpParam(row,p_idx,1,-Aineq(row,i));
                end

                % Check if parameter is used in Aeq --> add to beq
                rows = find(Aeq(:,i));
                for row=rows'
                    qcqpParams.beq(end+1) = newAdditiveQcqpParam(row,p_idx,1,-Aeq(row,i));
                end

                % Check bounds
                if lb(i) ~= -Inf || ub(i) ~= Inf
                    beep
                    warning('Y2F:parameterBounds','Bounds on parameters have no effect.')
                end

                % Quadratic constraints don't have to be checked

            elseif ~is(var, 'linear')

                % Nonlinear variables cannot be binary
                if binary
                    error('Nonlinear variables cannot be binary.');
                end

                % Is variable bilinear/quadratic?
                if sum(internalmodel.monomtable(i,:)) == 2
                    % Is variable the square of another variable?
                    if any(internalmodel.monomtable(i,:)==2)
                        v_idx = find(internalmodel.monomtable(i,:)==2,1);
                        assert(length(v_idx)==1);
                        v = internalmodel.used_variables(v_idx);
                        if any(paramVars == v)
                            error('Parameters can only be used affinely.')
                        end

                        % Remove this pseudo-variable
                        removeIdx = [removeIdx i];

                        % If it appears in the linear cost --> move it
                        if f(i) ~= 0
                            H(v_idx,v_idx) = H(v_idx,v_idx) + 2*f(i);
                        end

                        % Cannot appear in quadratic cost
                        if nnz(H(i,:)) > 0 || nnz(H(:,i)) > 0
                            error('Non-quadratic term appears in cost.')
                        end

                        % Cannot appear in equalities
                        if nnz(Aeq(:,i)) > 0
                            error('Quadratic equalities are not supported.')
                        end

                        % Check bounds
                        if lb(i) ~= -Inf || ub(i) ~= Inf
                            beep
                            warning('Y2F:quadraticTermBounds','Bounds on quadratic terms have no effect.')
                        end

                        % Check if variable is used in Aineq --> put in quad. ineq.
                        rows = find(Aineq(:,i));
                        for row=rows'
                            [k,quadIneq,Q,l,r] = findOrCreateQuadraticInequality(row,quadIneq,Q,l,r);
                            Q{k}(v_idx,v_idx) = Aineq(row,i);
                        end

                    % Is it bilinearly dependent!
                    % (e.g. if C_i is a parameter)
                    else
                        deps_idx = find(internalmodel.monomtable(i,:)==1);
                        assert(length(deps_idx) == 2);
                        deps = internalmodel.used_variables(deps_idx);

                        if all(ismember(deps, paramVars))
                            error('Parameters can only be used affinely.')
                        end

                        if ~any(ismember(deps, paramVars)) % no parameter
                            % Remove this pseudo-variable
                            removeIdx = [removeIdx i];

                            % If it appears in the linear cost --> move it
                            if f(i) ~= 0
                                H(deps_idx(1),deps_idx(2)) = H(deps_idx(1),deps_idx(2)) + f(i);
                                H(deps_idx(2),deps_idx(1)) = H(deps_idx(2),deps_idx(1)) + f(i);
                            end

                            % Cannot appear in quadratic cost
                            if nnz(H(i,:)) > 0 || nnz(H(:,i)) > 0
                                error('Non-quadratic term appears in cost.')
                            end

                            % Cannot appear in equalities
                            if nnz(Aeq(:,i)) > 0
                                error('Bilinear equalities are not supported.')
                            end

                            % Check bounds
                            if lb(i) ~= -Inf || ub(i) ~= Inf
                                beep
                                warning('Y2F:bilinearTermBounds','Bounds on bilinear terms have no effect.')
                            end

                            % Check if variable is used in Aineq --> put in quad. ineq.
                            rows = find(Aineq(:,i));
                            for row=rows'
                                [k,quadIneq,Q,l,r] = findOrCreateQuadraticInequality(row,quadIneq,Q,l,r);
                                Q{k}(deps_idx(1),deps_idx(2)) = 0.5*Aineq(row,i);
                                Q{k}(deps_idx(2),deps_idx(1)) = 0.5*Aineq(row,i);
                            end

                        else % a parameter is involed
                            if any(deps(1) == paramVars) && ...
                                    ~any(deps(2) == paramVars) && ...
                                    is(recover(deps(2)),'linear')
                                v = deps(2); % variable
                                v_idx = find(internalmodel.used_variables == v,1);
                                p = deps(1); % parameter
                                p_idx = find(paramVars == p,1);
                            elseif any(deps(2) == paramVars) && ...
                                    ~any(deps(1) == paramVars) && ...
                                    is(recover(deps(1)),'linear')
                                v = deps(1);
                                v_idx = find(internalmodel.used_variables == v,1);
                                p = deps(2);
                                p_idx = find(paramVars == p,1);
                            else % This shouldn't happen! Case is handled above
                                error('Internal error! Please contact support@embotech.ch')
                            end

                            % Variable is parameter --> remove it
                            removeIdx = [removeIdx i];

                            % Does bilinear combo influence cost?
                            if f(i) ~= 0
                                qcqpParams.f(end+1) = newAdditiveQcqpParam(v_idx,p_idx,1,f(i));
                            end

                            if any(H(:,i)) || any(H(i,:))
                                error('Parameters can only be used affinely.')
                            end

                            % Check if parameter is used in Aineq
                            rows = find(Aineq(:,i));
                            for row=rows'
                                qcqpParams.Aineq(end+1) = newAdditiveQcqpParam(sub2ind(size(Aineq),row,v_idx),p_idx,1,Aineq(row,i));
                            end

                            % Check if parameter is used in Aeq
                            rows = find(Aeq(:,i));
                            for row=rows'
                                qcqpParams.Aeq(end+1) = newAdditiveQcqpParam(sub2ind(size(Aeq),row,v_idx),p_idx,1,Aeq(row,i));
                            end

                            % Check bounds
                            if lb(i) ~= -Inf || ub(i) ~= Inf
                                beep
                                warning('Y2F:bilinearTermBounds','Bounds on bilinear terms have no effect.')
                            end
                        end
                    end
                elseif sum(internalmodel.monomtable(i,:))==3
                    % Decide which variable is the parameter
                    % At the end v1 and v2 are real variables, p is the parameter
                    temp_idx = find(internalmodel.monomtable(i,:));
                    if length(temp_idx) == 2 % p*x^2
                        v1_idx = find(internalmodel.monomtable(i,:) == 2,1); % variable
                        v1 = internalmodel.used_variables(v1_idx);
                        v2_idx = v1_idx;
                        v2 = v1;
                        p = internalmodel.used_variables(find(internalmodel.monomtable(i,:) == 1,1));
                        p_idx = find(paramVars == p,1);
                        if isempty(p_idx) || any(paramVars == v1)
                            error('One of the non-quadratic terms cannot be interpreted.')
                        end
                    elseif length(temp_idx) == 3 % p*x1*x2
                        temp_idx = find(internalmodel.monomtable(i,:) == 1);
                        temp_vars = internalmodel.used_variables(temp_idx); % variables & parameter
                        if ismember(temp_vars(1), paramVars) && ~any(ismember(temp_vars(2:3), paramVars))
                            v1_idx = temp_idx(2); % variables
                            v1 = temp_vars(2);
                            v2_idx = temp_idx(3);
                            v2 = temp_vars(3);
                            p = temp_vars(1);
                            p_idx = find(paramVars == p,1);
                        elseif ismember(temp_vars(2), paramVars) && ~any(ismember(temp_vars([1 3]), paramVars))
                            v1_idx = temp_idx(1); % variables
                            v1 = temp_vars(1);
                            v2_idx = temp_idx(3);
                            v2 = temp_vars(3);
                            p = temp_vars(2);
                            p_idx = find(paramVars == p,1);
                        elseif ismember(temp_vars(3), paramVars) && ~any(ismember(temp_vars(1:2), paramVars))
                            v1_idx = temp_idx(1); % variables
                            v1 = temp_vars(1);
                            v2_idx = temp_idx(2);
                            v2 = temp_vars(2);
                            p = temp_vars(3);
                            p_idx = find(paramVars == p,1);
                        else
                            error('One of the non-quadratic terms cannot be interpreted.')
                        end
                    else
                        error('One of the non-quadratic terms cannot be interpreted.')
                    end

                    % Variable is parameter --> remove it
                    removeIdx = [removeIdx i];

                    if f(i) ~= 0 % pseudo-variable affects cost --> param in H
                        qcqpParams.H(end+1) = newAdditiveQcqpParam(sub2ind(size(H),v1_idx,v2_idx),p_idx,1,f(i));
                        qcqpParams.H(end+1) = newAdditiveQcqpParam(sub2ind(size(H),v2_idx,v1_idx),p_idx,1,f(i));
                    end

                    % Cannot appear in quadratic cost
                    if nnz(H(i,:)) > 0 || nnz(H(:,i)) > 0
                        error('Non-quadratic term appears in cost.')
                    end

                    % Cannot appear in equalities
                    if nnz(Aeq(:,i)) > 0
                        error('Nonlinear equalities are not allowed.')
                    end

                    % Check bounds
                    if lb(i) ~= -Inf || ub(i) ~= Inf
                        beep
                        warning('Y2F:nonlinearTermBounds','Bounds on nonlinear terms have no effect.')
                    end

                    % Check inequalities
                    rows = find(Aineq(:,i));
                    for row=rows'
                        [k,quadIneq,Q,l,r] = findOrCreateQuadraticInequality(row,quadIneq,Q,l,r);
                        if v1_idx ~= v2_idx % make sure Q is symmetric
                            qcqpParams.Q(end+1) = newAdditiveQcqpParam(sub2ind(size(Q{k}),v1_idx,v2_idx),p_idx,k,0.5*Aineq(row,i));
                            qcqpParams.Q(end+1) = newAdditiveQcqpParam(sub2ind(size(Q{k}),v2_idx,v1_idx),p_idx,k,0.5*Aineq(row,i));
                        else
                            qcqpParams.Q(end+1) = newAdditiveQcqpParam(sub2ind(size(Q{k}),v1_idx,v1_idx),p_idx,k,Aineq(row,i));
                        end
                    end
                else
                    error('One of the non-quadratic terms cannot be interpreted.')
                end    
            end
        end

        % Before removing params, every real state can be recovered
        solverVars = internalmodel.used_variables;
        solverVars(removeIdx) = [];

        % Compute shifts of variables (to adjust parameters)
        shift = zeros(1,length(internalmodel.used_variables)); 
        for i=1:length(shift)
            shift(i) = nnz(removeIdx <= i);
        end

        % Remove parameter indices from equations
        H(removeIdx,:) = [];
        H(:,removeIdx) = [];
        f(removeIdx) = [];
        Aineq(:,removeIdx) = [];
        Aeq(:,removeIdx) = [];
        lb(removeIdx) = [];
        ub(removeIdx) = [];
        if ~isempty(Q)
            for k=1:numel(Q)
                Q{k}(removeIdx,:) = [];
                Q{k}(:,removeIdx) = [];
            end
            l(removeIdx,:) = [];
        end

        % Shift indices of parameters
        for i=1:numel(qcqpParams.H)
            [row,col] = ind2sub(size(H)+length(removeIdx),qcqpParams.H(i).maps2index);
            qcqpParams.H(i).maps2index = sub2ind(size(H),row-shift(row),col-shift(col));
        end
        for i=1:numel(qcqpParams.f)
            qcqpParams.f(i).maps2index = qcqpParams.f(i).maps2index - shift(qcqpParams.f(i).maps2index);
        end
        for i=1:numel(qcqpParams.Aineq)
            [row,col] = ind2sub(size(Aineq)+[0 length(removeIdx)],qcqpParams.Aineq(i).maps2index);
            qcqpParams.Aineq(i).maps2index = sub2ind(size(Aineq),row,col-shift(col));
        end
        for i=1:numel(qcqpParams.Aeq)
            [row,col] = ind2sub(size(Aeq)+[0 length(removeIdx)],qcqpParams.Aeq(i).maps2index);
            qcqpParams.Aeq(i).maps2index = sub2ind(size(Aeq),row,col-shift(col));
        end
        for i=1:numel(qcqpParams.lb)
            qcqpParams.lb(i).maps2index = qcqpParams.lb(i).maps2index - shift(qcqpParams.lb(i).maps2index);
        end
        for i=1:numel(qcqpParams.ub)
            qcqpParams.ub(i).maps2index = qcqpParams.ub(i).maps2index - shift(qcqpParams.ub(i).maps2index);
        end
        for i=1:numel(qcqpParams.Q)
            [row,col] = ind2sub(size(Q{qcqpParams.Q(i).maps2mat})+length(removeIdx),qcqpParams.Q(i).maps2index);
            qcqpParams.Q(i).maps2index = sub2ind(size(Q{qcqpParams.Q(i).maps2mat}),row,col-shift(col));
        end
        for i=1:numel(qcqpParams.l)
            qcqpParams.l(i).maps2index = qcqpParams.l(i).maps2index - shift(qcqpParams.l(i).maps2index);
        end

        % Shift & sort binary variables
        qcqpParams.bidx = sort(qcqpParams.bidx - shift(qcqpParams.bidx));

        % Finish converting linear inequalities to quadratic ones
        for i=1:size(quadIneq,2)
            row = quadIneq(1,i);
            k = quadIneq(2,i);

            r(k) = r(k) + bineq(row);
            l(:,k) = l(:,k) + Aineq(row,:)';

            % Convert bineq params
            relevantParams = findRelevantParams(row, 1, size(bineq), qcqpParams.bineq);
            for j=fliplr(relevantParams) % sort descending
                qcqpParams.r(end+1) = newAdditiveQcqpParam(k,qcqpParams.bineq(j).maps2origparam, ...
                                        1, qcqpParams.bineq(j).factor);
                qcqpParams.bineq(j) = [];
            end

            % Convert Aineq params
            relevantParams = findRelevantParams(row, 1:size(Aineq,2), size(Aineq), qcqpParams.Aineq);
            for j=fliplr(relevantParams) % sort descending
                [~,col] = ind2sub(size(Aineq), qcqpParams.Aineq(j).maps2index);
                qcqpParams.l(end+1) = newAdditiveQcqpParam(col,qcqpParams.Aineq(j).maps2origparam, ...
                                        k, qcqpParams.Aineq(j).factor);
                qcqpParams.Aineq(j) = [];
            end
        end
        % Compute shift in rows for linear inequalites
        shift = zeros(1,length(bineq));
        for i=1:length(shift)
            shift(i) = nnz(quadIneq(1,:) <= i);
        end
        % Delete linear inequalities
        Aineq(quadIneq(1,:),:) = [];
        bineq(quadIneq(1,:)) = [];
        % Fix parameters
        for i=1:numel(qcqpParams.Aineq)
            [row,col] = ind2sub(size(Aineq)+[length(quadIneq(1,:)) 0],qcqpParams.Aineq(i).maps2index);
            qcqpParams.Aineq(i).maps2index = sub2ind(size(Aineq),row-shift(row),col);
        end
        for i=1:numel(qcqpParams.bineq)
            qcqpParams.bineq(i).maps2index = qcqpParams.bineq(i).maps2index - shift(qcqpParams.bineq(i).maps2index);
        end

        % Convert inequalities to bounds - Round 2
        % This is only necessary if we have parameters that affect bounds
        bounds_idx = [];
        for i=1:size(Aineq,1)
            vars = find(Aineq(i,:)); % variables used in inequality
            if length(vars) == 1 && is(recover(solverVars(vars)),'linear')
                relevantParams = findRelevantParams(i, 1:size(Aineq,2), size(Aineq), qcqpParams.Aineq);
                if isempty(relevantParams) % No parameters affecting LHS of inequality
                    if Aineq(i,vars) > 0 && ub(vars) == Inf % No bounds so far
                        % we can make a bound out of this
                        bounds_idx(end+1) = i;
                        ub(vars) = bineq(i)/Aineq(i,vars);

                        % Convert bineq params
                        relevantParams = findRelevantParams(i, 1, size(bineq), qcqpParams.bineq);
                        for j=fliplr(relevantParams) % sort descending
                            qcqpParams.ub(end+1) = newAdditiveQcqpParam(vars,qcqpParams.bineq(j).maps2origparam, ...
                                                    1, qcqpParams.bineq(j).factor/Aineq(i,vars));
                            qcqpParams.bineq(j) = [];
                        end
                    elseif Aineq(i,vars) < 0 && lb(vars) == -Inf % No bounds so far
                        % we can make a bound out of this
                        bounds_idx(end+1) = i;
                        lb(vars) = bineq(i)/Aineq(i,vars);

                        % Convert bineq params
                        relevantParams = findRelevantParams(i, 1, size(bineq), qcqpParams.bineq);
                        for j=fliplr(relevantParams) % sort descending
                            qcqpParams.lb(end+1) = newAdditiveQcqpParam(vars,qcqpParams.bineq(j).maps2origparam, ...
                                                    1, qcqpParams.bineq(j).factor/Aineq(i,vars));
                            qcqpParams.bineq(j) = [];
                        end
                    end
                end
            end
        end
        % Compute shift in rows for linear inequalites
        shift = zeros(1,length(bineq));
        for i=1:length(shift)
            shift(i) = nnz(bounds_idx <= i);
        end
        % Delete rows
        Aineq(bounds_idx,:) = [];
        bineq(bounds_idx) = [];
        % Fix parameters
        for i=1:numel(qcqpParams.Aineq)
            [row,col] = ind2sub(size(Aineq)+[length(bounds_idx) 0],qcqpParams.Aineq(i).maps2index);
            qcqpParams.Aineq(i).maps2index = sub2ind(size(Aineq),row-shift(row),col);
        end
        for i=1:numel(qcqpParams.bineq)
            qcqpParams.bineq(i).maps2index = qcqpParams.bineq(i).maps2index - shift(qcqpParams.bineq(i).maps2index);
        end
    end

    function [k,quadIneq,Q,l,r] = findOrCreateQuadraticInequality(rowIdx,quadIneq,Q,l,r)
    % Helper function to find/create quadratic inequalities
    %   Input:
    %       rowIdx      index of linear inequality (row in Aineq)
    %       quadIneq    data structure containing information about quadratic
    %                   inequalities
    %                   first row:     id of linear inequality
    %                   second row:    id of quad inequality
    %       Q,l,r       quadratic inequalities
    %
    %   Output:
    %       k           index of quadratic inequality (as in Q{k} = ...)
    %       updated quadIneq, Q, l, and r
        id = find(quadIneq(1,:) == rowIdx);
        if ~isempty(id)
            assert(length(id) == 1);
            k = quadIneq(2,id);
        else
            Q{end+1} = spalloc(size(H,1),size(H,2),0);%zeros(size(H));
            l(:,end+1) = spalloc(size(H,1),1,0); %zeros(size(H,1),1);
            r(end+1) = 0;
            k = length(r);
            quadIneq(:,end+1) = [rowIdx;k];
        end
    end

    function outputFORCES = buildOutput()
    % Helper function that builds the output struct required for the FORCES
    % solver(s), an outputMap that allows to recover the wantend output 
    % values from the solver output, an outputParamTable that allows the
    % usage of parameters in outputs
    
        outputFORCES = {};
        % we need to know which output to get from which solver
        sys.outputMap = zeros(3,0); % 1st row: variable type (1=decision variable,2=parameter)
                                    % 2nd row: index of solver/index of parameter value
                                    % 3rd row: index of output/index of element inside value matrix

        o = ones(numel(stages),1); % counter variable outputs
        p = 1; % counter parameters
        k = 1; % counter total number of outputs
        for i=1:numel(solverOutputs)
            outputVars = getvariables(solverOutputs{i});
            sys.outputBase{i} = full(getbase(solverOutputs{i}));
            sys.outputSize{i} = size(solverOutputs{i});
            for j=1:length(outputVars)
                idx = find(solverVars == outputVars(j),1);
                if length(idx) == 1
                    [stage, state, component] = findVariableIndex(graphComponents,idx);
                    outputFORCES{component}(o(component)) = newOutput(sprintf('o_%u',o(component)), stage, state);
                    sys.outputMap(:,end+1) = [1; component; o(component)];
                    o(component) = o(component) + 1;
                else
                    idx = find(paramVars == outputVars(j),1);
                    if length(idx) == 1
                        sys.outputMap(:,end+1) = [2; yalmipParamMap(:,idx)];
                        p = p+1;
                    else
                        error('Output is not valid. Only linear combinations of optimization variables and parameters are allowed.')
                    end
                end
                k = k + 1;
            end
        end
        sys.lengthOutput = k-1;
    end

    function checkQcqpForInfeasibility()
    % Helper function that checks if the QCQP might be infeasible and warns
    % the user
        
        % Check if cost matrix is positive semi-definite
        if isempty(qcqpParams.H) && (~issymmetric(H) || ~all(eig(H) >= -1e-7))
            error('Hessian is not positive semi-definite.')
        end

        if ~issymmetric(H_temp)
            beep
            warning('Y2F:nonsymmetricHessian','Hessian is not symmetric.')
        end

        % Check if bounds make sense
        if isempty(qcqpParams.lb) && isempty(qcqpParams.ub) && any(lb > ub)
            error('Bounds are infeasible.')
        end

        % Check if any quadratic constraint matrix is indefinite
        for k=1:numel(Q_temp)
            if ~issymmetric(Q_temp{k})
                beep
                warning('Y2F:nonsymmetricQuadraticConstraint', ...
                    'One of the quadratic constraint matrices is not symmetric.')
            end

            if isempty(qcqpParams.Q) && any(eig(Q{k}) < -1e-7) && any(eig(Q{k}) > -1e-7)
                beep
                warning('Y2F:indefiniteQuadraticConstraint', ...
                    'One of the quadratic constraint matrices is indefinite.')
            end
        end
        
        % Create warning if binary variables are used
        if ~isempty(qcqpParams.bidx)
            while 1
                beep
                % MATLAB displays anything between [\8 ... ]\8 in orange font
                in = input(['[' 8 'YALMIP created a mixed integer problem. This could slow down the generated solver.\nDo you want to continue [y/n]? ]' 8],'s');
                if strcmpi(in,'y')
                    break
                elseif strcmpi(in,'n')
                    error('Please reformulate your problem or contact support at support@embotech.com.');
                end
            end
        end
    end

end

