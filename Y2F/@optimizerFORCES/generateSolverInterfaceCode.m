function success = generateSolverInterfaceCode( self )
%GENERATESOLVERMANAGMENTCODE generates MEX C code that will prepare the
%parameters for the FORCES solver. It also assembles the correct outputs.

% Create file
solverName = self.codeoptions{1}.name;
if ~isdir(solverName)
    error('Solver ''%s'' has not been generated!', solverName)
end
fileID = fopen([solverName '/interface/' solverName '_optimizerFORCES_interface_mex.c'],'w');

% Write standard comment
fprintf(fileID, '/*\n This is an interface for %s that ',solverName);
fprintf(fileID, 'is used by optimizerFORCES to call the solver\n');
fprintf(fileID, '*/ \n\n');

% Includes
fprintf(fileID, '#include "mex.h"\n');
fprintf(fileID, '#include "math.h"\n');
for k=1:self.numSolvers
    fprintf(fileID, '#include "../../%s/include/%s.h"\n',self.codeoptions{k}.name,self.codeoptions{k}.name);
end
fprintf(fileID, '#include <stdio.h>\n\n');

% Copy functions stolen from FORCES
fprintf(fileID, '/* copy functions */\n');
fprintf(fileID, 'void copyCArrayToM(double *src, double *dest, int dim) {\n');
fprintf(fileID, '\twhile (dim--) {\n');
fprintf(fileID, '\t\t*dest++ = (double)*src++;\n');
fprintf(fileID, '\t}\n');
fprintf(fileID, '}\n');
fprintf(fileID, 'void copyMArrayToC(double *src, double *dest, int dim) {\n');
fprintf(fileID, '\twhile (dim--) {\n');
fprintf(fileID, '\t\t*dest++ = (double) (*src++) ;\n');
fprintf(fileID, '\t}\n');
fprintf(fileID, '}\n\n');

% Arguments for solver(s)
fprintf(fileID, '/* Some memory for mex-function */\n');
for k=1:self.numSolvers
    if ~self.solverIsBinary(k) % no binary variables
        fprintf(fileID, '%s_params params_%u;\n',self.codeoptions{k}.name,k);
        fprintf(fileID, '%s_output output_%u;\n',self.codeoptions{k}.name,k);
        fprintf(fileID, '%s_info info_%u;\n\n',self.codeoptions{k}.name,k);
    else
        fprintf(fileID, '%s_binaryparams params_%u;\n',self.codeoptions{k}.name,k);
        fprintf(fileID, '%s_binaryoutput output_%u;\n',self.codeoptions{k}.name,k);
        fprintf(fileID, '%s_info info_%u;\n\n',self.codeoptions{k}.name,k);
    end
end

% Start of MEX function
fprintf(fileID, '/* THE mex-function */\n');
fprintf(fileID, 'void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {\n');
fprintf(fileID, '\t/* file pointer for printing */\n');
fprintf(fileID, '\tFILE *fp = NULL;\n\n');
    
fprintf(fileID, '\t/* define variables */\n');
fprintf(fileID, '\tconst mxArray *param_values = prhs[0];\n');
fprintf(fileID, '\tmxArray *outvar;\n');
fprintf(fileID, '\tint i;\n');
fprintf(fileID, '\tint exitflag_1');
for i=2:self.numSolvers
    fprintf(fileID, ', exitflag_%u',i);
end
fprintf(fileID, ';\n');
fprintf(fileID, '\tconst char *fname;\n');

% Solver info fields
fprintf(fileID, '\tconst char *infofields[16] = { "it", "it2opt", "res_eq", "res_ineq",  "pobj",  "dobj",  "dgap", "rdgap",  "mu",  "mu_aff",  "sigma",  "lsit_aff",  "lsit_cc",  "step_aff",   "step_cc",  "solvetime"};\n\n');

% Check number of in- and outputs
fprintf(fileID, '\t/* Check for proper number of arguments */\n');
fprintf(fileID, '\tif (nrhs != 1) {\n');
fprintf(fileID, '\t\tmexErrMsgTxt("This function requires exactly 1 input: parameter value cell array.");\n');
fprintf(fileID, '\t}\n');
fprintf(fileID, '\tif (nlhs > 3) {\n');
fprintf(fileID, '\t\tmexErrMsgTxt("This function returns at most 3 outputs.");\n');
fprintf(fileID, '\t}\n');

% Check type of input
fprintf(fileID, '\t/* Check whether params is actually a cell array */\n');
fprintf(fileID, '\tif( !mxIsCell(param_values) ) {\n');
fprintf(fileID, '\t\tmexErrMsgTxt("param_values must be a cell array.");\n');
fprintf(fileID, '\t}\n\n');

% Check length of input
fprintf(fileID, '\t/* Check whether params is actually a cell array */\n');
fprintf(fileID, '\tif( mxGetNumberOfElements(param_values) != %u ) {\n',self.numParams);
fprintf(fileID, '\t\tmexErrMsgTxt("param_values must have %u elements.");\n',self.numParams);
fprintf(fileID, '\t}\n\n');

% Check sizes of parameter values and load their values
fprintf(fileID, '\t/* Load param values into C arrays and check their size */\n');
fprintf(fileID, '\tmxArray *par;\n\n');
for i=1:self.numParams
    % Get cell element
    fprintf(fileID, '\tpar = mxGetCell(param_values,%u);\n',i-1);
    
    % Cell element not found
    fprintf(fileID, '\tif( par == NULL ) {\n');
    fprintf(fileID, '\t\tmexErrMsgTxt("Parameter #%u not found");\n',i);
    fprintf(fileID, '\t}\n');
    
    % Parameter value is not double
    fprintf(fileID, '\tif( !mxIsDouble(par) ) {\n');
    fprintf(fileID, '\t\tmexErrMsgTxt("Parameter #%u must be a double.");\n',i);
    fprintf(fileID, '\t}\n');
    
    % Check parameter value size
    fprintf(fileID, '\tif( mxGetM(par) != %u || mxGetN(par) != %u ) {\n',self.paramSizes(i,1),self.paramSizes(i,2));
    fprintf(fileID, '\t\tmexErrMsgTxt("Parameter #%u must be of size [%u x %u]");\n',i,self.paramSizes(i,1),self.paramSizes(i,2));
    fprintf(fileID, '\t}\n');
    
    % Get parameter value
    fprintf(fileID, '\tconst double *param_value_%u = mxGetPr(par);\n\n',i);
end

% We need managment code for every solver
for k=1:self.numSolvers
    solverName = self.codeoptions{k}.name;
    
	% Set FORCES parameter values
    fprintf(fileID, '\t/*Assigning parameter values of solver #%u*/\n',k);
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
                    fprintf(fileID, '\tparams_%u.%s[%u] = %.15g;\n',k,fields{i},j-1,problem.(fields{i})(j));
                else
                    if abs(problem.(fields{i})(j)) > 1e-15 % only print standard value if it's not 0
                        fprintf(fileID, '\tparams_%u.%s[%u] = %.15g\n',k,fields{i},j-1,problem.(fields{i})(j)); % print std value
                        
                        % Add additive parameters
                        for p=ps
                            factor = paramMap(2,p);
                            valueMatrix = paramMap(3,p);
                            valueIndex = paramMap(4,p);
                            fprintf(fileID, ' + %.15g * param_value_%u[%u]',factor,valueMatrix,valueIndex-1);
                        end
                        fprintf(fileID, ';\n');
                    else
                        % print first additive param
                        factor = paramMap(2,ps(1));
                        valueMatrix = paramMap(3,ps(1));
                        valueIndex = paramMap(4,ps(1));
                        fprintf(fileID, '\tparams_%u.%s[%u] = %.15g * param_value_%u[%u]',k,fields{i},j-1,factor,valueMatrix,valueIndex-1);
                        % Add other additive parameters
                        for p=ps(2:end)
                            factor = paramMap(2,p);
                            valueMatrix = paramMap(3,p);
                            valueIndex = paramMap(4,p);
                            fprintf(fileID, ' + %.15g * param_value_%u[%u]',factor,valueMatrix,valueIndex-1);
                        end
                        fprintf(fileID, ';\n');
                    end
                end
            end
            fprintf(fileID, '\n');
        end
    else % we need the value of the fake parameter
        for i=1:numel(self.standardParamValues{k})
            fprintf(fileID, '\tparams_%u.p[%u] = %.15g;\n',k,i-1,self.standardParamValues{k}(i));
        end
    end
    
    % If the user wanted output, we need to store it
    fprintf(fileID, '\t#if %s_SET_PRINTLEVEL > 0\n',solverName);
    fprintf(fileID, '\t\t/* Prepare file for printfs */\n');
    fprintf(fileID, '\t\t/*fp = freopen("stdout_temp","w+",stdout);*/\n');
    fprintf(fileID, '\t\tfp = fopen("stdout_temp","w+");\n');
    fprintf(fileID, '\t\tif( fp == NULL ) {\n');
    fprintf(fileID, '\t\t\tmexErrMsgTxt("freopen of stdout did not work.");\n');
    fprintf(fileID, '\t\t}\n');
    fprintf(fileID, '\t\trewind(fp);\n');
	fprintf(fileID, '\t#endif\n\n');

	fprintf(fileID, '\t/* call solver #%u */\n',k);
	fprintf(fileID, '\texitflag_%u = %s_solve(&params_%u, &output_%u, &info_%u, fp );\n\n',k,solverName,k,k,k);

	fprintf(fileID, '\t/* close stdout */\n');
	fprintf(fileID, '\t/* fclose(fp); */\n\n');
	
    % Print output in console
	fprintf(fileID, '\t#if %s_SET_PRINTLEVEL > 0\n',solverName);
    fprintf(fileID, '\t\t/* Read contents of printfs printed to file */\n');
    fprintf(fileID, '\t\trewind(fp);\n');
    fprintf(fileID, '\t\twhile( (i = fgetc(fp)) != EOF ) {\n');
    fprintf(fileID, '\t\t\tmexPrintf("%%c",i);\n');
    fprintf(fileID, '\t\t}\n');
    fprintf(fileID, '\t\tfclose(fp);\n');
	fprintf(fileID, '\t#endif\n\n');
end
    
% Put outputs together
fprintf(fileID, '\t/*Build outputs*/\n');
mapOffset = 0; % needed for outputMap
if self.outputIsCell % multiple outputs possible
    % Create MATLAB output cell array
    fprintf(fileID, '\tplhs[0] = mxCreateCellMatrix(1, %u);\n',numel(self.outputBase));
    for i=1:numel(self.outputBase) % every output has a base
        base = self.outputBase{k};
        outSize = self.outputSize{k};
        fprintf(fileID, '\toutvar = mxCreateDoubleMatrix(%u, %u, mxREAL);\n',outSize(1),outSize(2));
        fprintf(fileID, '\tdouble *output_value_%u = mxGetPr(outvar);\n',i);
        %fprintf(fileID, '\tdouble output_value_%u[%u];\n',i,size(base,1));
        for j=1:size(base,1) % elements (decision var or parameters) needed for this output
            fprintf(fileID, '\toutput_value_%u[%u] = %.15g',i,j-1,base(j,1));
            idx = find(base(j,2:end)); % index of necessary elements inside outputMap (plus offset)
            for k=idx
                if self.outputMap(1,k+mapOffset) == 1
                    fprintf(fileID, ' + %.15g * output_%u.o_%u[0]',base(j,k+1),self.outputMap(2,k+mapOffset),self.outputMap(3,k+mapOffset));
                elseif self.outputMap(1,k+mapOffset) == 2
                    fprintf(fileID, ' + %.15g * param_value_%u[%u]',base(j,k+1),self.outputMap(2,k+mapOffset),self.outputMap(3,k+mapOffset)-1);
                else
                    error('Unknown output variable type');
                end
            end

            fprintf(fileID, ';\n');
        end
        %fprintf(fileID, '\tcopyCArrayToM(output_value_%u, mxGetPr(outvar), %u);\n',i,prod(outSize));
        fprintf(fileID, '\tmxSetCell(plhs[0], %u, outvar);\n\n',i-1);
        mapOffset = mapOffset + size(base,1);
    end
else % only one output
    % Create MATLAB output matrix
    base = self.outputBase{1};
    outSize = self.outputSize{1};
    fprintf(fileID, '\tplhs[0] = mxCreateDoubleMatrix(%u, %u, mxREAL);\n',outSize(1),outSize(2));
    fprintf(fileID, '\tdouble *output_value = mxGetPr(plhs[0]);\n');
    %fprintf(fileID, '\tdouble output_value[%u];\n',size(base,1));
    for j=1:size(base,1) % elements (decision var or parameters) needed for this output
        fprintf(fileID, '\toutput_value[%u] = %.15g',j-1,base(j,1));
        idx = find(base(j,2:end)); % index of necessary elements inside outputMap (plus offset)
        for k=idx
            if self.outputMap(1,k+mapOffset) == 1
                fprintf(fileID, ' + %.15g * output_%u.o_%u[0]',base(j,k+1),self.outputMap(2,k+mapOffset),self.outputMap(3,k+mapOffset));
            elseif self.outputMap(1,k+mapOffset) == 2
                fprintf(fileID, ' + %.15g * param_value_%u[%u]',base(j,k+1),self.outputMap(2,k+mapOffset),self.outputMap(3,k+mapOffset)-1);
            else
                error('Unknown output variable type');
            end
        end
        
        fprintf(fileID, ';\n\n');
    end
    %fprintf(fileID, '\tcopyCArrayToM(output_value, mxGetPr(plhs[0]), %u);\n\n',prod(outSize));
end

% add other output arguments (exitflags and info)
if self.numSolvers > 1
    fprintf(fileID, '\t/* copy exitflags */\n');
    fprintf(fileID, '\tif( nlhs > 1 ) {\n');
    fprintf(fileID, '\t\tplhs[1] = mxCreateDoubleMatrix(1, %u, mxREAL);\n',self.numSolvers);
    for i=1:self.numSolvers
        fprintf(fileID, '\t\t*(mxGetPr(plhs[1])+%u) = (double)exitflag_%u;\n',i-1,i);
    end
    fprintf(fileID, '\t}\n\n');

    fprintf(fileID, '\t/* copy info structs */\n');
    fprintf(fileID, '\tif( nlhs > 2 ) {\n');
    fprintf(fileID, '\t\tplhs[2] = mxCreateCellMatrix(1, %u);\n',self.numSolvers);
    fprintf(fileID, '\t\tmxArray *temp_info;\n\n');

    for i=1:self.numSolvers
        fprintf(fileID, '\t\t/* info for solver #%u */\n',i);
        fprintf(fileID, '\t\ttemp_info = mxCreateStructMatrix(1, 1, 16, infofields);\n\n');

        fprintf(fileID, '\t\t/* iterations */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_%u.it;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "it", outvar);\n\n');

        fprintf(fileID, '\t\t/* iterations to optimality (branch and bound) */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_%u.it2opt;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "it2opt", outvar);\n\n');

        fprintf(fileID, '\t\t/* res_eq */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.res_eq;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "res_eq", outvar);\n\n');

        fprintf(fileID, '\t\t/* res_ineq */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.res_ineq;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "res_ineq", outvar);\n\n');

        fprintf(fileID, '\t\t/* pobj */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.pobj;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "pobj", outvar);\n\n');

        fprintf(fileID, '\t\t/* dobj */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.dobj;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "dobj", outvar);\n\n');

        fprintf(fileID, '\t\t/* dgap */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.dgap;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "dgap", outvar);\n\n');

        fprintf(fileID, '\t\t/* rdgap */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.rdgap;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "rdgap", outvar);\n\n');

        fprintf(fileID, '\t\t/* mu */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.mu;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "mu", outvar);\n\n');

        fprintf(fileID, '\t\t/* mu_aff */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.mu_aff;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "mu_aff", outvar);\n\n');

        fprintf(fileID, '\t\t/* sigma */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.sigma;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "sigma", outvar);\n\n');

        fprintf(fileID, '\t\t/* lsit_aff */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_%u.lsit_aff;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "lsit_aff", outvar);\n\n');

        fprintf(fileID, '\t\t/* lsit_cc */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_%u.lsit_cc;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "lsit_cc", outvar);\n\n');

        fprintf(fileID, '\t\t/* step_aff */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.step_aff;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "step_aff", outvar);\n\n');

        fprintf(fileID, '\t\t/* step_cc */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.step_cc;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "step_cc", outvar);\n\n');

        fprintf(fileID, '\t\t/* solver time */\n');
        fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
        fprintf(fileID, '\t\t*mxGetPr(outvar) = info_%u.solvetime;\n',i);
        fprintf(fileID, '\t\tmxSetField(temp_info, 0, "solvetime", outvar);\n\n');

        fprintf(fileID, '\t\tmxSetCell(plhs[2],%u,temp_info);\n\n',i-1);
    end
    fprintf(fileID, '\t}\n');
else % only one solver --> no cell arrays necessary
    fprintf(fileID, '\t/* copy exitflag */\n');
    fprintf(fileID, '\tif( nlhs > 1 ) {\n');
    fprintf(fileID, '\t\tplhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*(mxGetPr(plhs[1])) = (double)exitflag_1;\n');
    fprintf(fileID, '\t}\n\n');

    fprintf(fileID, '\t/* copy info struct */\n');
    fprintf(fileID, '\tif( nlhs > 2 ) {\n');
    fprintf(fileID, '\t\tplhs[2] = mxCreateStructMatrix(1, 1, 16, infofields);\n\n');

    fprintf(fileID, '\t\t/* iterations */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_1.it;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "it", outvar);\n\n');

    fprintf(fileID, '\t\t/* iterations to optimality (branch and bound) */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_1.it2opt;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "it2opt", outvar);\n\n');

    fprintf(fileID, '\t\t/* res_eq */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.res_eq;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "res_eq", outvar);\n\n');

    fprintf(fileID, '\t\t/* res_ineq */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.res_ineq;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "res_ineq", outvar);\n\n');

    fprintf(fileID, '\t\t/* pobj */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.pobj;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "pobj", outvar);\n\n');

    fprintf(fileID, '\t\t/* dobj */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.dobj;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "dobj", outvar);\n\n');

    fprintf(fileID, '\t\t/* dgap */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.dgap;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "dgap", outvar);\n\n');

    fprintf(fileID, '\t\t/* rdgap */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.rdgap;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "rdgap", outvar);\n\n');

    fprintf(fileID, '\t\t/* mu */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.mu;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "mu", outvar);\n\n');

    fprintf(fileID, '\t\t/* mu_aff */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.mu_aff;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "mu_aff", outvar);\n\n');

    fprintf(fileID, '\t\t/* sigma */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.sigma;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "sigma", outvar);\n\n');

    fprintf(fileID, '\t\t/* lsit_aff */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_1.lsit_aff;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "lsit_aff", outvar);\n\n');

    fprintf(fileID, '\t\t/* lsit_cc */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = (double)info_1.lsit_cc;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "lsit_cc", outvar);\n\n');

    fprintf(fileID, '\t\t/* step_aff */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.step_aff;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "step_aff", outvar);\n\n');

    fprintf(fileID, '\t\t/* step_cc */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.step_cc;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "step_cc", outvar);\n\n');

    fprintf(fileID, '\t\t/* solver time */\n');
    fprintf(fileID, '\t\toutvar = mxCreateDoubleMatrix(1, 1, mxREAL);\n');
    fprintf(fileID, '\t\t*mxGetPr(outvar) = info_1.solvetime;\n');
    fprintf(fileID, '\t\tmxSetField(plhs[2], 0, "solvetime", outvar);\n');
    
    fprintf(fileID, '\t}\n');
end

fprintf(fileID, '}'); % end of mex-function
    
% Don't forget to close MEX file
fclose(fileID);

success = 1;

end

