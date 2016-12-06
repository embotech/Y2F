function compileSolverInterfaceCode( self )
%COMPILESOLVERINTERFACECODE Compiles the MEX code generated by GENERATESOLVERINTERFACECODE

solverName = self.default_codeoptions.name;
cName = [solverName '/interface/' solverName];
mexName = [solverName '/interface/' solverName '_mex'];
outputName = ['"' solverName '"'];
    


% copy the O-files of all solvers into /interface
% we'll delete them later, but this makes compilation easier
for i=1:self.numSolvers
    if( ispc )
        copyfile(sprintf('%s/obj/%s.obj',self.codeoptions{i}.name,self.codeoptions{i}.name), sprintf('%s/interface',solverName), 'f');
    else % mac or linux
        copyfile(sprintf('%s/obj/%s.o',self.codeoptions{i}.name,self.codeoptions{i}.name), sprintf('%s/interface',solverName), 'f');
    end
end

if exist( [cName '.c'], 'file' ) && exist( [mexName '.c'], 'file' )
    mex('-c','-O','-outdir',[solverName '/interface'],[cName '.c'])
    mex('-c','-O','-outdir',[solverName '/interface'],[mexName '.c'])
    if( ispc )
        mex([solverName '/interface/*.obj'], '-output', outputName) 
        delete([solverName '/interface/*.obj']);
    elseif( ismac )
        mex([solverName '/interface/*.o'], '-output', outputName) 
        delete([solverName '/interface/*.o']);
    else % we're on a linux system
        mex([solverName '/interface/*.o'], '-output', outputName,'-lrt') 
        delete([solverName '/interface/*.o']);
    end
else
    fprintf('Could not find source file. This file is meant to be used for building from source code.');
end

end

