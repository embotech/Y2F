function compileSolverInterfaceCode( self )
%COMPILESOLVERINTERFACECODE Compiles the MEX code generated by
%GENERATECINTERFACECODE and GENERATEMEXINTERFACECODE.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech AG, Zurich, Switzerland, 2013-2020.

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

% copy the object/library files of all solvers into /interface
% we'll delete them later, but this makes compilation easier
if( ~ispc )
    missingObjs = false;
    missingLibs = false;
    useLibraryFiles = true;
    % check whether libraries and object files exist
    for i=1:self.numSolvers
        objName = sprintf('%s/obj/%s.o',solverName,self.codeoptions{i}.name); 
        libName = sprintf('%s/lib/lib%s.a',solverName,self.codeoptions{i}.name);
        if(~exist(objName, 'file'))
            missingObjs = true;
        end
        if(~exist(libName, 'file'))
            missingLibs = true;
        end
    end
    if(missingLibs && missingObjs)
        error('Cannot find object files or libraries to compile');
    end
    if(~missingLibs)
        useLibraryFiles = true;
        for i=1:self.numSolvers
            copyfile(sprintf('%s/lib/lib%s.a',solverName,self.codeoptions{i}.name), sprintf('%s/interface',solverName), 'f');
        end
    else
        useLibraryFiles = false;
        for i=1:self.numSolvers
            copyfile(sprintf('%s/obj/%s.o',solverName,self.codeoptions{i}.name), sprintf('%s/interface',solverName), 'f');
        end
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
        intelLibsDir = [clientPath,filesep,'libs_Intel'];
        if( exist( intelLibsDir, 'dir' ) && ~isMinGW)
            % If subdirectory for specific architecture exisst, we need to
            % use that
            archDir = [intelLibsDir, filesep, 'win', arch];
            if( exist(archDir, 'dir') )
                intelLibsDir = archDir;
            end
            intelLibsDirFlag = ['-L', intelLibsDir];
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
                iphlpapi_path = ['-L"', matlabroot, '\sys\lcc64\lcc64\lib64"'];
                mex([solverName '/interface/' solverName '.obj'], ...
                    [solverName '/interface/' solverName '_mex.obj'], ...
                    '-output', outputName, ...
                    [solverName '/lib/', libs{:}], ...
                    '-lIPHLPAPI.lib', iphlpapi_path, '-largeArrayDims', '-silent');
            end
        else
            % it seems that we have been compiling with VS only,
            % so we do not add the Intel libs and use only object files
            
            % we create a cell array for mex interface, main solver and all internal solvers
            linkFiles = cell(1,self.numSolvers + 2);
            linkFiles{1} = [solverName, '/interface/', solverName, '.obj'];
            linkFiles{2} = [solverName, '/interface/', solverName, '_mex.obj'];
            % Add all object files in obj folder
            for i=1:self.numSolvers
                objFile = [solverName, '/obj/', self.codeoptions{i}.name, '.obj'];
                if exist(objFile, 'file')
                    linkFiles{2 + i} = objFile;
                else
                    error(['Object file ', objFile, ' not found. Compilation aborted');
                end
            end

            % Compile MEX interface
            mex(linkFiles{:}, '-lIPHLPAPI.lib', legacyLibs, ...
                '-output', [outputName(2:end-1),'.',mexext], '-largeArrayDims', '-silent');
        
            % Delete unnecessary object files
            delete([solverName '/interface/*.obj']);
        end
    else % macOS or linux system
        
        % we create a cell array for mex interface, main solver and all internal solvers
        linkFiles = cell(1,self.numSolvers + 2);
        linkFiles{1} = [solverName, '/interface/', solverName, '.o'];
        linkFiles{2} = [solverName, '/interface/', solverName, '_mex.o'];
        if(useLibraryFiles)
            % Add all library files in lib folder
            for i=1:self.numSolvers
                libFile = [solverName, '/lib/lib', self.codeoptions{i}.name, '.a'];
                linkFiles{2 + i} = libFile;
            end
        else
            % Add all object files in obj folder
            for i=1:self.numSolvers
                objFile = [solverName, '/obj/', self.codeoptions{i}.name, '.o'];
                linkFiles{2 + i} = objFile;
            end
        end    
        
        % Compile MEX interface
        if( ismac ) % macOS system
            mex(linkFiles{:}, '-output', outputName, '-largeArrayDims', '-silent');
        else % linux system
            mex(linkFiles{:}, '-output', outputName, '-lrt', '-largeArrayDims', '-silent');
        end
        
        % Delete unnecessary library/object files
        if(useLibraryFiles)
            delete([solverName '/interface/lib*.a']);
        end
        delete([solverName '/interface/*.o']);
    end
else
    fprintf('Could not find source file. This file is meant to be used for building from source code.');
end

end

