function compileSolverInterfaceCode( self )
%COMPILESOLVERINTERFACECODE Compiles the MEX code generated by
%GENERATECINTERFACECODE and GENERATEMEXINTERFACECODE.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech AG, Zurich, Switzerland, 2013-2019.

solverName = self.default_codeoptions.name;
cName = [solverName '/interface/' solverName];
mexName = [solverName '/interface/' solverName '_mex'];
outputName = ['"' solverName '"'];
    
% move the (necessary) files of all solvers to the new directory and delete
% the folders of the "internal" solvers
for i=1:self.numSolvers
    % include
    dir2move = sprintf('%s/include',self.codeoptions{i}.name);
    if( exist(dir2move,'dir') )
        copyfile(dir2move, sprintf('%s/include',solverName), 'f');
    end
    % lib
    dir2move = sprintf('%s/lib',self.codeoptions{i}.name);
    if( exist(dir2move,'dir') )
        copyfile(dir2move, sprintf('%s/lib',solverName), 'f');
    end
    % obj
    dir2move = sprintf('%s/obj',self.codeoptions{i}.name);
    if( exist(dir2move,'dir') )
        copyfile(dir2move, sprintf('%s/obj',solverName), 'f');
    end
    % src
    dir2move = sprintf('%s/src',self.codeoptions{i}.name);
    if exist(dir2move,'dir')
        copyfile(dir2move, sprintf('%s/src',solverName), 'f');
    end
    % obj_target
    dir2move = sprintf('%s/obj_target',self.codeoptions{i}.name);
    if exist(dir2move,'dir')
        copyfile(dir2move, sprintf('%s/obj_target',solverName), 'f');
    end
    % lib_target
    dir2move = sprintf('%s/lib_target',self.codeoptions{i}.name);
    if exist(dir2move,'dir')
        copyfile(dir2move, sprintf('%s/lib_target',solverName), 'f');
    end
    
    % Delete files
    rmdir(self.codeoptions{i}.name, 's');
    delete([self.codeoptions{i}.name '*']);
end

try 
    thisCompiler = mex.getCompilerConfigurations('C','Selected');
    mexcomp.name = thisCompiler(1).Name;
    mexcomp.ver = thisCompiler(1).Version;
    mexcomp.vendor = thisCompiler(1).Manufacturer;
catch
    mexcomp = [];
end

if(isstruct(mexcomp) && isfield(mexcomp, 'name') && strncmpi(mexcomp.name, 'MinGW', 5))
    isMinGW = true;
else
    isMinGW = false;
end

% copy the O-files of all solvers into /interface
% we'll delete them later, but this makes compilation easier
for i=1:self.numSolvers
    if( ~ispc )
        copyfile(sprintf('%s/obj/%s.o',solverName,self.codeoptions{i}.name), sprintf('%s/interface',solverName), 'f');
    end
end

% final MEX build
if exist( [cName '.c'], 'file' ) && exist( [mexName '.c'], 'file' )
    if isfield(self.default_codeoptions, 'optlevel') && self.default_codeoptions.optlevel > 0
        optFlags = '-O'; % enable code optimization
    else
        optFlags = '-g'; % disable optimization by enabling debug flags
    end
    
    mex('-c',optFlags,'-largeArrayDims', '-silent','-outdir',[solverName '/interface'],[cName '.c']) % compiles C interface
    mex('-c',optFlags,'-largeArrayDims', '-silent','-outdir',[solverName '/interface'],[mexName '.c']) % compiles MEX interface
    
    if( ispc ) % PC - we need additional libraries
        
        % Create a list of internal solver libraries
        if( exist([solverName,filesep,'lib'],'dir') )
            libs = cell(1,self.numSolvers);
            if(isMinGW)
                for i=1:self.numSolvers
                    lib = dir([solverName,filesep,'lib/lib',self.codeoptions{i}.name,'*.a']);
                    if length(lib)>1
                        % fix for new server which produces shared and static
                        % libraries, so we have more than 1 library found above
                        lib = dir([solverName,filesep,'lib/lib',self.codeoptions{i}.name,'.a']);
                    end
                    libs{i} = [lib.name];
                end
            else
                for i=1:self.numSolvers
                    lib = dir([solverName,filesep,'lib/',self.codeoptions{i}.name,'*.lib']);
                    if length(lib)>1
                        % fix for new server which produces shared and static
                        % libraries, so we have more than 1 library found above
                        lib = dir([solverName,filesep,'lib/',self.codeoptions{i}.name,'_static.lib']);
                    end
                    libs{i} = ['-l' lib.name(1:end-4)];
                end
            end
        end
        
        % figure out if we are on a 32 or 64 bit system
        ext = mexext;
        arch = ext(end-1:end); % mex extension ends with 32 or 64
        
        % figure our whether we need additional libraries for Intel
        clientPath = fileparts(which('generateCode'));
        intelLibsDir = [clientPath,filesep,'libs_intel'];
        if( exist( intelLibsDir, 'dir' ) && ~isMinGW)
            % If subdirectory for specific architecture exisst, we need to
            % use that
            archDir = [intelLibsDir, filesep, 'win', arch];
            if( exist(archDir, 'dir') )
                intelLibsDir = archDir;
            end
            intelLibsDirFlag = ['-L', intelLibsDir];
            addpath(intelLibsDir); savepath;
        else
            intelLibsDirFlag = '';
        end
        
        % Figure out whether we need legacy libraries for Visual Studio        
        try
            thisCompiler = mex.getCompilerConfigurations('C','Selected');
            mexcomp.name = thisCompiler(1).Name;
            mexcomp.ver = thisCompiler(1).Version;
            mexcomp.vendor = thisCompiler(1).Manufacturer;
        catch
            mexcomp = [];
        end
        if( ~isempty(mexcomp) ...
           && ~isempty(strfind(mexcomp.vendor,'Microsoft')) ...
           && str2double(mexcomp.ver) >= 14.0 )
           legacyLibs = '-llegacy_stdio_definitions';
        else
           legacyLibs = '';
        end
        
        % Call mex compiler
        if( exist([solverName,filesep,'lib'],'dir') )
            if(~isMinGW)
                mex([solverName '/interface/' solverName '.obj'], ...
                    [solverName '/interface/' solverName '_mex.obj'], ...
                    '-output', outputName, ...
                    ['-L' solverName '/lib'], libs{:}, intelLibsDirFlag, ...
                    '-ldecimal', '-lirc', '-lmmt', '-lsvml_dispmt', ...
                    legacyLibs, '-lIPHLPAPI.lib', '-largeArrayDims', '-silent');
            else
                mex([solverName '/interface/' solverName '.obj'], ...
                    [solverName '/interface/' solverName '_mex.obj'], ...
                    '-output', outputName, ...
                    [solverName '/lib/', libs{:}], ...
                    '-lIPHLPAPI.lib', '-largeArrayDims', '-silent');
            end
        else
            % it seems that we have been compiling with VS only,
            % so we do not add the Intel libs and use only object files
            
            % Find all object files in obj folder
            objFiles = dir([solverName filesep 'obj']); % struct with name and folder in different fields
            objFiles = objFiles(~cellfun(@isempty, regexp({objFiles.name}, '\.obj$', 'start', 'once')));
            objFiles = arrayfun(@(x) [solverName filesep 'obj' filesep x.name], ...
                                objFiles, 'UniformOutput', false); % cell array with paths

            % Compile MEX interface
            mex([solverName, '/interface/' solverName '.obj'], ...
                [solverName, '/interface/' solverName '_mex.obj'], ...
                objFiles{:}, '-lIPHLPAPI.lib', legacyLibs, ...
                '-output', [outputName(2:end-1),'.',mexext], '-largeArrayDims', '-silent');
        
            % Delete unnecessary object files
            delete([solverName '/interface/*.obj']);
        end
    elseif( ismac ) % macOS
        
        % Find all object files in interface folder
        objFiles = dir([solverName filesep 'interface']); % struct with name and folder in different fields
        objFiles = objFiles(~cellfun(@isempty, regexp({objFiles.name}, '\.o$', 'start', 'once')));
        objFiles = arrayfun(@(x) [solverName filesep 'interface' filesep x.name], ...
                            objFiles, 'UniformOutput', false); % cell array with paths
        
        % Compile MEX interface
        mex(objFiles{:}, '-output', outputName, '-largeArrayDims', '-silent')
        
        % Delete unnecessary object files
        delete([solverName '/interface/*.o']);
    else % we're on a linux system
        
        % Find all object files in interface folder
        objFiles = dir([solverName filesep 'interface']); % struct with name and folder in different fields
        objFiles = objFiles(~cellfun(@isempty, regexp({objFiles.name}, '\.o$', 'start', 'once')));
        objFiles = arrayfun(@(x) [solverName filesep 'interface' filesep x.name], ...
                            objFiles, 'UniformOutput', false); % cell array with paths
        
        % Compile MEX interface
        mex(objFiles{:}, '-output', outputName, '-lrt', '-largeArrayDims', '-silent') 
        
        % Delete unnecessary object files
        delete([solverName '/interface/*.o']);
    end
else
    fprintf('Could not find source file. This file is meant to be used for building from source code.');
end

end

