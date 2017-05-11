function success = generateCInterfaceCode( self )
%GENERATECINTERFACECODE generates C code that will prepare the user-defined
%parameters for the FORCES solver. It also assembles the correct outputs.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech GmbH, Zurich, Switzerland, 2013-2016.

success = 0;

% Get solver name from option
solverName = self.default_codeoptions.name;

% Check if FORCES solver has been generated
if ~isdir(self.codeoptions{1}.name) && ~isdir(solverName)
    error('Solver ''%s'' has not been generated!', solverName)
end

% Make directories
if ~isdir(solverName)
    mkdir(solverName)
end
if ~isdir([solverName '/interface'])
    mkdir([solverName '/interface'])
end
if ~isdir([solverName '/include'])
    mkdir([solverName '/include'])
end
if ~isdir([solverName '/solvers'])
    mkdir([solverName '/solvers'])
end

% % Move "internal" solver(s) to new directory ("hide" them from the user)
% for i=1:self.numSolvers
%     copyfile(self.codeoptions{i}.name, [solverName '/solvers/' self.codeoptions{i}.name]);
% end

% Remove "internal" solver(s)
% for i=1:self.numSolvers
%     FORCEScleanup(self.codeoptions{i}.name, 'all');
% end


% Create files
hFileID = fopen([solverName '/include/' solverName '.h'],'w');
cFileID = fopen([solverName '/interface/' solverName '.c'],'w');

% Create h file
fprintf(hFileID, '/*\nHeader file containing definitions for C interface of %s,\n', solverName);
fprintf(hFileID, ' a fast costumized optimization solver.\n*/\n\n');

fprintf(hFileID, '#include <stdio.h>\n\n');

fprintf(hFileID, '#ifndef __%s_H__\n',solverName);
fprintf(hFileID, '#define __%s_H__\n\n',solverName);

% Visual Studio 2015 Compatibility
fprintf(hFileID, '/* For Visual Studio 2015 Compatibility */\n');
fprintf(hFileID, '#if _MSC_VER == 1900\n');
fprintf(hFileID, 'FILE * __cdecl __iob_func(void);\n');
fprintf(hFileID, '#endif\n');

fprintf(hFileID, '/* DATA TYPE ------------------------------------------------------------*/\n');
fprintf(hFileID, 'typedef double %s_FLOAT;\n\n',solverName);

%fprintf(hFileID, 'typedef double %sINTERFACE_FLOAT;\n\n',solverName);

fprintf(hFileID, '/* SOLVER SETTINGS ------------------------------------------------------*/\n');
fprintf(hFileID, '/* print level */\n');
fprintf(hFileID, '#ifndef %s_SET_PRINTLEVEL\n', solverName);
fprintf(hFileID, '#define %s_SET_PRINTLEVEL    (%u)\n', solverName, self.default_codeoptions.printlevel);
fprintf(hFileID, '#endif\n\n');

fprintf(hFileID, '/* PARAMETERS -----------------------------------------------------------*/\n');
fprintf(hFileID, '/* fill this with data before calling the solver! */\n');
fprintf(hFileID, 'typedef struct %s_params\n',solverName);
fprintf(hFileID, '{\n');

for i=1:self.numParams
    if self.paramSizes(i,1) == 1 && self.paramSizes(i,2) == 1 % scalar parameter
        fprintf(hFileID, '\t/* scalar */\n');
    elseif self.paramSizes(i,1) == 1 % row vector
        fprintf(hFileID, '\t/* row vector of length %u */\n',self.paramSizes(i,2));
    elseif self.paramSizes(i,2) == 1 % column vector
        fprintf(hFileID, '\t/* column vector of length %u */\n',self.paramSizes(i,1));
    else
        fprintf(hFileID, '\t/* matrix of size [%u x %u] (column major format) */\n',self.paramSizes(i,1),self.paramSizes(i,2));
    end
    fprintf(hFileID, '\t%s_FLOAT %s[%u];\n\n',solverName,self.paramNames{i},prod(self.paramSizes(i,:)));
end

fprintf(hFileID, '} %s_params;\n\n\n',solverName);


fprintf(hFileID, '/* OUTPUTS --------------------------------------------------------------*/\n');
fprintf(hFileID, '/* the desired variables are put here by the solver */\n');
fprintf(hFileID, 'typedef struct %s_output\n',solverName);
fprintf(hFileID, '{\n');

for i=1:numel(self.outputSize)
    if self.outputSize{i}(1) == 1 && self.outputSize{i}(2) == 1 % scalar parameter
        fprintf(hFileID, '\t/* scalar */\n');
    elseif self.outputSize{i}(1) == 1 % row vector
        fprintf(hFileID, '\t/* row vector of length %u */\n',self.outputSize{i}(2));
    elseif self.outputSize{i}(2) == 1 % column vector
        fprintf(hFileID, '\t/* column vector of length %u */\n',self.outputSize{i}(1));
    else
        fprintf(hFileID, '\t/* matrix of size [%u x %u] (column major format) */\n',self.outputSize{i}(1),self.outputSize{i}(2));
    end
    fprintf(hFileID, '\t%s_FLOAT %s[%u];\n\n',solverName,self.outputNames{i},prod(self.outputSize{i}));
end

fprintf(hFileID, '} %s_output;\n\n\n',solverName);

% Solver info (contains arrays if multiple solvers are used)
if self.numSolvers == 1
    fprintf(hFileID, '/* SOLVER INFO ----------------------------------------------------------*/\n');
    fprintf(hFileID, '/* diagnostic data from last interior point step */\n');
    fprintf(hFileID, 'typedef struct %s_info\n',solverName);
    fprintf(hFileID, '{\n');
    
    fprintf(hFileID, '\t/* iteration number */\n');
    fprintf(hFileID, '\tint it;\n\n');

	fprintf(hFileID, '\t/* number of iterations needed to optimality (branch-and-bound) */\n');
	fprintf(hFileID, '\tint it2opt;\n\n');
	
    fprintf(hFileID, '\t/* inf-norm of equality constraint residuals */\n');
    fprintf(hFileID, '\t%s_FLOAT res_eq;\n\n',solverName);
	
    if isfield(self.default_codeoptions, 'solvemethod') && strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % extra field for ADMM
        fprintf(hFileID, '\t/* inf-norm of dual residual */\n');
        fprintf(hFileID, '\t%s_FLOAT res_dual;\n\n',solverName);
    else % not for ADMM
        fprintf(hFileID, '\t/* inf-norm of inequality constraint residuals */\n');
        fprintf(hFileID, '\t%s_FLOAT res_ineq;\n\n',solverName);
    end

    fprintf(hFileID, '\t/* primal objective */\n');
    fprintf(hFileID, '\t%s_FLOAT pobj;\n\n',solverName);
	
    fprintf(hFileID, '\t/* dual objective */\n');
    fprintf(hFileID, '\t%s_FLOAT dobj;\n\n',solverName);

    fprintf(hFileID, '\t/* duality gap := pobj - dobj */\n');
    fprintf(hFileID, '\t%s_FLOAT dgap;\n\n',solverName);
	
    fprintf(hFileID, '\t/* relative duality gap := |dgap / pobj | */\n');
    fprintf(hFileID, '\t%s_FLOAT rdgap;\n\n',solverName);

    if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
        fprintf(hFileID, '\t/* duality measure */\n');
        fprintf(hFileID, '\t%s_FLOAT mu;\n\n',solverName);
        
        fprintf(hFileID, '\t/* duality measure (after affine step) */\n');
        fprintf(hFileID, '\t%s_FLOAT mu_aff;\n\n',solverName);
        
        fprintf(hFileID, '\t/* centering parameter */\n');
        fprintf(hFileID, '\t%s_FLOAT sigma;\n\n',solverName);
        
        fprintf(hFileID, '\t/* number of backtracking line search steps (affine direction) */\n');
        fprintf(hFileID, '\tint lsit_aff;\n\n');
        
        fprintf(hFileID, '\t/* number of backtracking line search steps (combined direction) */\n');
        fprintf(hFileID, '\tint lsit_cc;\n\n');
        
        fprintf(hFileID, '\t/* step size (affine direction) */\n');
        fprintf(hFileID, '\t%s_FLOAT step_aff;\n\n',solverName);
        
        fprintf(hFileID, '\t/* step size (combined direction) */\n');
        fprintf(hFileID, '\t%s_FLOAT step_cc;\n\n',solverName);
    end

	fprintf(hFileID, '\t/* solvertime */\n');
	fprintf(hFileID, '\t%s_FLOAT solvetime;\n\n',solverName);

    fprintf(hFileID, '} %s_info;\n\n\n',solverName);
else
    fprintf(hFileID, '/* SOLVER INFO ----------------------------------------------------------*/\n');
    fprintf(hFileID, '/* diagnostic data from last interior point step for every solver */\n');
    fprintf(hFileID, '/* (in total %u solvers are used) */\n',self.numSolvers);
    fprintf(hFileID, 'typedef struct %s_info\n',solverName);
    fprintf(hFileID, '{\n');
    
    fprintf(hFileID, '\t/* iteration number */\n');
    fprintf(hFileID, '\tint it[%u];\n\n',self.numSolvers);

	fprintf(hFileID, '\t/* number of iterations needed to optimality (branch-and-bound) */\n');
	fprintf(hFileID, '\tint it2opt[%u];\n\n',self.numSolvers);
	
    fprintf(hFileID, '\t/* inf-norm of equality constraint residuals */\n');
    fprintf(hFileID, '\t%s_FLOAT res_eq[%u];\n\n',solverName,self.numSolvers);
    
    if isfield(self.default_codeoptions, 'solvemethod') && strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % extra field for ADMM
        fprintf(hFileID, '\t/* inf-norm of inequality constraint residuals */\n');
        fprintf(hFileID, '\t%s_FLOAT res_dual[%u];\n\n',solverName,self.numSolvers);
    else % not for ADMM
        fprintf(hFileID, '\t/* inf-norm of inequality constraint residuals */\n');
        fprintf(hFileID, '\t%s_FLOAT res_ineq[%u];\n\n',solverName,self.numSolvers);
    end

    fprintf(hFileID, '\t/* primal objective */\n');
    fprintf(hFileID, '\t%s_FLOAT pobj[%u];\n\n',solverName,self.numSolvers);
	
    fprintf(hFileID, '\t/* dual objective */\n');
    fprintf(hFileID, '\t%s_FLOAT dobj[%u];\n\n',solverName,self.numSolvers);

    fprintf(hFileID, '\t/* duality gap := pobj - dobj */\n');
    fprintf(hFileID, '\t%s_FLOAT dgap[%u];\n\n',solverName,self.numSolvers);
	
    fprintf(hFileID, '\t/* relative duality gap := |dgap / pobj | */\n');
    fprintf(hFileID, '\t%s_FLOAT rdgap[%u];\n\n',solverName,self.numSolvers);

    if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
        fprintf(hFileID, '\t/* duality measure */\n');
        fprintf(hFileID, '\t%s_FLOAT mu[%u];\n\n',solverName,self.numSolvers);

        fprintf(hFileID, '\t/* duality measure (after affine step) */\n');
        fprintf(hFileID, '\t%s_FLOAT mu_aff[%u];\n\n',solverName,self.numSolvers);

        fprintf(hFileID, '\t/* centering parameter */\n');
        fprintf(hFileID, '\t%s_FLOAT sigma[%u];\n\n',solverName,self.numSolvers);

        fprintf(hFileID, '\t/* number of backtracking line search steps (affine direction) */\n');
        fprintf(hFileID, '\tint lsit_aff[%u];\n\n',self.numSolvers);

        fprintf(hFileID, '\t/* number of backtracking line search steps (combined direction) */\n');
        fprintf(hFileID, '\tint lsit_cc[%u];\n\n',self.numSolvers);

        fprintf(hFileID, '\t/* step size (affine direction) */\n');
        fprintf(hFileID, '\t%s_FLOAT step_aff[%u];\n\n',solverName,self.numSolvers);

        fprintf(hFileID, '\t/* step size (combined direction) */\n');
        fprintf(hFileID, '\t%s_FLOAT step_cc[%u];\n\n',solverName,self.numSolvers);
    end

	fprintf(hFileID, '\t/* solvertime */\n');
	fprintf(hFileID, '\t%s_FLOAT solvetime[%u];\n\n',solverName,self.numSolvers);

    fprintf(hFileID, '} %s_info;\n\n\n',solverName);
end

% Solver function
fprintf(hFileID, '/* SOLVER FUNCTION DEFINITION -------------------------------------------*/\n');
if self.numSolvers == 1
    fprintf(hFileID, '/* examine exitflag before using the result! */\n');
    fprintf(hFileID, 'int %s_solve(%s_params* params, %s_output* output, %s_info* info, FILE* fs);\n\n', solverName, solverName, solverName, solverName);
else
    fprintf(hFileID, '/* examine all of the %u exitflags before using the result! */\n', self.numSolvers);
    fprintf(hFileID, 'int* %s_solve(%s_params* params, %s_output* output, %s_info* info, FILE* fs);\n\n', solverName, solverName, solverName, solverName);
end

fprintf(hFileID, '#endif\n');

% Close h-file
fclose(hFileID);

% Write standard comment
fprintf(cFileID, '/*\n This is an interface for %s that ',solverName);
fprintf(cFileID, 'can be used to call the solver generated by FORCES Pro\n');
fprintf(cFileID, '*/ \n\n');

fprintf(cFileID, '#include "../include/%s.h"\n',solverName);
for k=1:self.numSolvers
    fprintf(cFileID, '#include "../include/%s.h"\n',self.codeoptions{k}.name);
end
fprintf(cFileID, '#include <stdio.h>\n\n');

% Visual Studio 2015 Compatibility
fprintf(cFileID, '/* For Visual Studio 2015 Compatibility */\n');
fprintf(cFileID, '#if _MSC_VER == 1900\n');
fprintf(cFileID, 'FILE _iob[3];\n');
fprintf(cFileID, 'FILE * __cdecl __iob_func(void)\n');
fprintf(cFileID, '{\n');
fprintf(cFileID, '    _iob[0] = *stdin;\n');
fprintf(cFileID, '    _iob[1] = *stdout;\n');
fprintf(cFileID, '    _iob[2] = *stderr;\n');
fprintf(cFileID, '    return _iob;\n');
fprintf(cFileID, '}\n');
fprintf(cFileID, '#endif\n\n');

% Arguments for solver(s)
fprintf(cFileID, '/* Some memory */\n');
for k=1:self.numSolvers
    if ~self.solverIsBinary(k) % no binary variables
        fprintf(cFileID, '%s_params params_%u;\n',self.codeoptions{k}.name,k);
        fprintf(cFileID, '%s_output output_%u;\n',self.codeoptions{k}.name,k);
        fprintf(cFileID, '%s_info info_%u;\n\n',self.codeoptions{k}.name,k);
    else
        fprintf(cFileID, '%s_binaryparams params_%u;\n',self.codeoptions{k}.name,k);
        fprintf(cFileID, '%s_binaryoutput output_%u;\n',self.codeoptions{k}.name,k);
        fprintf(cFileID, '%s_info info_%u;\n\n',self.codeoptions{k}.name,k);
    end
end
if self.numSolvers == 1
    fprintf(cFileID, 'int exitflag;\n');
else
    fprintf(cFileID, 'int exitflag[%u];\n',self.numSolvers);
end

% Start of interface function
if self.numSolvers == 1
    fprintf(cFileID, 'int %s_solve(%s_params* params, %s_output* output, %s_info* info, FILE* fs) {\n', solverName, solverName, solverName, solverName);
else
    fprintf(cFileID, 'int* %s_solve(%s_params* params, %s_output* output, %s_info* info, FILE* fs) {\n', solverName, solverName, solverName, solverName);
end
    
fprintf(cFileID, '\t/* define variables */\n');
fprintf(cFileID, '\tint i;\n');

% We need managment code for every solver
for k=1:self.numSolvers
    fprintf(cFileID, '\t/* SOLVER %u --------------------------------------------------------*/\n',k);
    
    fprintf(cFileID, '\t/*Assigning parameter values of solver #%u*/\n',k);
    solverName = self.codeoptions{k}.name;
    
	% Set FORCES parameter values
    fprintf(cFileID, '\t/*Assigning parameter values of solver #%u*/\n',k);
    if self.solverHasParams(k) % otherwise there is fake one
        problem = self.standardParamValues{k};
        fields = fieldnames(problem);
        % Go through all parameters (fields contains p_1, p_2, ...)
        for i=1:numel(fields)
            paramMap = self.forcesParamMap{k}.(fields{i});
            % Go through all elements of this parameter
            for j=1:numel(problem.(fields{i}))
                ps = find(paramMap(1,:) == j); % all relevant param map entries
                if isempty(ps)
                    % Load standard value
                    fprintf(cFileID, '\tparams_%u.%s[%u] = %.15g;\n',k,fields{i},j-1,problem.(fields{i})(j));
                else
                    if abs(problem.(fields{i})(j)) > 1e-15 % only print standard value if it's not 0
                        fprintf(cFileID, '\tparams_%u.%s[%u] = %.15g',k,fields{i},j-1,problem.(fields{i})(j)); % print std value
                        
                        % Add additive parameters
                        for p=ps
                            factor = paramMap(2,p);
                            valueMatrix = paramMap(3,p);
                            valueIndex = paramMap(4,p);
                            fprintf(cFileID, ' + %.15g * params->%s[%u]',factor,self.paramNames{valueMatrix},valueIndex-1);
                        end
                        fprintf(cFileID, ';\n');
                    else
                        % print first additive param
                        factor = paramMap(2,ps(1));
                        valueMatrix = paramMap(3,ps(1));
                        valueIndex = paramMap(4,ps(1));
                        fprintf(cFileID, '\tparams_%u.%s[%u] = %.15g * params->%s[%u]',k,fields{i},j-1,factor,self.paramNames{valueMatrix},valueIndex-1);
                        % Add other additive parameters
                        for p=ps(2:end)
                            factor = paramMap(2,p);
                            valueMatrix = paramMap(3,p);
                            valueIndex = paramMap(4,p);
                            fprintf(cFileID, ' + %.15g * params->%s[%u]',factor,self.paramNames{valueMatrix},valueIndex-1);
                        end
                        fprintf(cFileID, ';\n');
                    end
                end
            end
            fprintf(cFileID, '\n');
        end
    else % we need the value of the fake parameter
        for i=1:numel(self.standardParamValues{k})
            fprintf(cFileID, '\tparams_%u.p[%u] = %.15g;\n',k,i-1,self.standardParamValues{k}(i));
        end
    end

	fprintf(cFileID, '\t/* call solver #%u */\n',k);
    if self.numSolvers == 1
        fprintf(cFileID, '\texitflag = %s_solve(&params_%u, &output_%u, &info_%u, fs );\n\n',solverName,k,k,k);
    else
        fprintf(cFileID, '\texitflag[%u] = %s_solve(&params_%u, &output_%u, &info_%u, fs );\n\n',k-1,solverName,k,k,k);
    end
    
    % Diagnostics
    if self.numSolvers == 1
        fprintf(cFileID, '\t/* iterations */\n');
        fprintf(cFileID, '\tinfo->it = info_%u.it;\n\n',k);

        fprintf(cFileID, '\t/* iterations to optimality (branch and bound) */\n');
        fprintf(cFileID, '\tinfo->it2opt = info_%u.it2opt;\n\n',k);

        fprintf(cFileID, '\t/* res_eq */\n');
        fprintf(cFileID, '\tinfo->res_eq = info_%u.res_eq;\n\n',k);

        if isfield(self.default_codeoptions, 'solvemethod') && strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % extra field for ADMM
            fprintf(cFileID, '\t/* res_dual */\n');
            fprintf(cFileID, '\tinfo->res_dual = info_%u.res_dual;\n\n',k);            
        else % not for ADMM
            fprintf(cFileID, '\t/* res_ineq */\n');
            fprintf(cFileID, '\tinfo->res_ineq = info_%u.res_ineq;\n\n',k);
        end

        fprintf(cFileID, '\t/* pobj */\n');
        fprintf(cFileID, '\tinfo->pobj = info_%u.pobj;\n\n',k);

        fprintf(cFileID, '\t/* dobj */\n');
        fprintf(cFileID, '\tinfo->dobj = info_%u.dobj;\n\n',k);

        fprintf(cFileID, '\t/* dgap */\n');
        fprintf(cFileID, '\tinfo->dgap = info_%u.dgap;\n\n',k);

        fprintf(cFileID, '\t/* rdgap */\n');
        fprintf(cFileID, '\tinfo->rdgap = info_%u.rdgap;\n\n',k);
        
        if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
            fprintf(cFileID, '\t/* mu */\n');
            fprintf(cFileID, '\tinfo->mu = info_%u.mu;\n\n',k);

            fprintf(cFileID, '\t/* mu_aff */\n');
            fprintf(cFileID, '\tinfo->mu_aff = info_%u.mu_aff;\n\n',k);

            fprintf(cFileID, '\t/* sigma */\n');
            fprintf(cFileID, '\tinfo->sigma = info_%u.sigma;\n\n',k);

            fprintf(cFileID, '\t/* lsit_aff */\n\n');
            fprintf(cFileID, '\tinfo->lsit_aff = info_%u.lsit_aff;\n\n',k);

            fprintf(cFileID, '\t/* lsit_cc */\n');
            fprintf(cFileID, '\tinfo->lsit_cc = info_%u.lsit_cc;\n\n',k);

            fprintf(cFileID, '\t/* step_aff */\n');
            fprintf(cFileID, '\tinfo->step_aff = info_%u.step_aff;\n\n',k);

            fprintf(cFileID, '\t/* step_cc */\n');
            fprintf(cFileID, '\tinfo->step_cc = info_%u.step_cc;\n\n',k);
        end

        fprintf(cFileID, '\t/* solver time */\n');
        fprintf(cFileID, '\tinfo->solvetime = info_%u.solvetime;\n\n\n',k);
    else 
        fprintf(cFileID, '\t/* iterations */\n');
        fprintf(cFileID, '\tinfo->it[%u] = info_%u.it;\n\n',k-1,k);

        fprintf(cFileID, '\t/* iterations to optimality (branch and bound) */\n');
        fprintf(cFileID, '\tinfo->it2opt[%u] = info_%u.it2opt;\n\n',k-1,k);

        fprintf(cFileID, '\t/* res_eq */\n');
        fprintf(cFileID, '\tinfo->res_eq[%u] = info_%u.res_eq;\n\n',k-1,k);

        
        if isfield(self.default_codeoptions, 'solvemethod') && strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % extra field for ADMM
            fprintf(cFileID, '\t/* res_ineq */\n');
            fprintf(cFileID, '\tinfo->res_dual[%u] = info_%u.res_dual;\n\n',k-1,k);
        else % not for ADMM
            fprintf(cFileID, '\t/* res_ineq */\n');
            fprintf(cFileID, '\tinfo->res_ineq[%u] = info_%u.res_ineq;\n\n',k-1,k);
        end

        fprintf(cFileID, '\t/* pobj */\n');
        fprintf(cFileID, '\tinfo->pobj[%u] = info_%u.pobj;\n\n',k-1,k);

        fprintf(cFileID, '\t/* dobj */\n');
        fprintf(cFileID, '\tinfo->dobj[%u] = info_%u.dobj;\n\n',k-1,k);

        fprintf(cFileID, '\t/* dgap */\n');
        fprintf(cFileID, '\tinfo->dgap[%u] = info_%u.dgap;\n\n',k-1,k);

        fprintf(cFileID, '\t/* rdgap */\n');
        fprintf(cFileID, '\tinfo->rdgap[%u] = info_%u.rdgap;\n\n',k-1,k);

        if ~isfield(self.default_codeoptions, 'solvemethod') || ~strcmpi(self.default_codeoptions.solvemethod, 'ADMM') % not for ADMM
            fprintf(cFileID, '\t/* mu */\n');
            fprintf(cFileID, '\tinfo->mu[%u] = info_%u.mu;\n\n',k-1,k);

            fprintf(cFileID, '\t/* mu_aff */\n');
            fprintf(cFileID, '\tinfo->mu_aff[%u] = info_%u.mu_aff;\n\n',k-1,k);

            fprintf(cFileID, '\t/* sigma */\n');
            fprintf(cFileID, '\tinfo->sigma[%u] = info_%u.sigma;\n\n',k-1,k);

            fprintf(cFileID, '\t/* lsit_aff */\n');
            fprintf(cFileID, '\tinfo->lsit_aff[%u] = info_%u.lsit_aff;\n\n',k-1,k);

            fprintf(cFileID, '\t/* lsit_cc */\n');
            fprintf(cFileID, '\tinfo->lsit_cc[%u] = info_%u.lsit_cc;\n\n',k-1,k);

            fprintf(cFileID, '\t/* step_aff */\n');
            fprintf(cFileID, '\tinfo->step_aff[%u] = info_%u.step_aff;\n\n',k-1,k);

            fprintf(cFileID, '\t/* step_cc */\n');
            fprintf(cFileID, '\tinfo->step_cc[%u] = info_%u.step_cc;\n\n',k-1,k);
        end

        fprintf(cFileID, '\t/* solver time */\n');
        fprintf(cFileID, '\tinfo->solvetime[%u] = info_%u.solvetime;\n\n\n',k-1,k);
    end
end
    
% Put outputs together
fprintf(cFileID, '\t/* OUTPUTS -----------------------------------------------------------*/\n');
fprintf(cFileID, '\t/*Build outputs*/\n');
mapOffset = 0; % needed for outputMap
for i=1:numel(self.outputBase) % every output has a base
    base = self.outputBase{i};
    for j=1:size(base,1) % elements (decision var or parameters) needed for this output
        fprintf(cFileID, '\toutput->%s[%u] = %.15g',self.outputNames{i},j-1,base(j,1));
        idx = find(base(j,2:end)); % index of necessary elements inside outputMap (plus offset)
        for k=idx
            if self.outputMap(1,k+mapOffset) == 1
                fprintf(cFileID, ' + %.15g * output_%u.o_%u[0]',base(j,k+1),self.outputMap(2,k+mapOffset),self.outputMap(3,k+mapOffset));
            elseif self.outputMap(1,k+mapOffset) == 2
                fprintf(cFileID, ' + %.15g * params->%s[%u]',base(j,k+1),self.paramNames{self.outputMap(2,k+mapOffset)},self.outputMap(3,k+mapOffset)-1);
            else
                error('Unknown output variable type');
            end
        end

        fprintf(cFileID, ';\n\n\n');
    end
    mapOffset = mapOffset + size(base,1);
end


fprintf(cFileID, '\treturn exitflag;\n');
fprintf(cFileID, '}'); % end of mex-function
    
% Don't forget to close C file
fclose(cFileID);

success = 1;

end

