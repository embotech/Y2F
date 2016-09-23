function result = checkIfGraphIsConnected( G )
%CHECKIFGRAPHISCONNECTED result = 1 if G is connected graph, 0 otherwise
% G has the format described in EMPTYGRAPH

if G.n <= 1
    result = 1;
    return
end

nextVertices = 1;
verticesNotVisited = 2:G.n;

while ~isempty(verticesNotVisited) && ~isempty(nextVertices)
    v = nextVertices(1);
    idx = intersect(find(G.adjMatrix(v,:)), verticesNotVisited);
    verticesNotVisited = setdiff(verticesNotVisited, idx);
    nextVertices = [nextVertices(2:end) idx];
end

result = isempty(verticesNotVisited);

end

