function adj = generateAdjMatrixFromEqualities( Aeq )
%GENERATEADJMATRIXFROMEQUALITIES Generates adjacency matrix from equality
%constraints. If two variables appear in the same equality, the
%corresponding vertices are connected.

n = size(Aeq,2);
adj = double(Aeq ~= 0); % temp result
adj = (adj'*adj) > 0;
adj(1:n+1:end) = 0; % diagonal has to be 0

end

