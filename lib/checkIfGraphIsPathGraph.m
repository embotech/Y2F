function result = checkIfGraphIsPathGraph( G )
%CHECKGRAPHISPATHGRAPH Checks if a graph is a path graph
%   result = 1 if given graph is a path graph, 0 otherwise
%   G has the format described in EMPTYGRAPH

result = 0;

if ~checkIfGraphIsConnected(G)
    return
end

% Simple cases
if G.n <= 1
    result = 1;
    return
end

% Do we have the right degrees? We need two 1s and the rest 2s
deg = sum(G.adjMatrix);
% Find the 1s
ones = find(deg == 1);
if length(ones) == 2
    deg(ones) = [];
    if all(deg == 2)
        result = 1;
    end
end

end

