function G = subgraph( G, idx )
%SUBGRAPH Returns the induced subgraph of G containing the vertices with
% indices in idx

G.vertices = G.vertices(idx);
G.adjMatrix = G.adjMatrix(idx, idx);
G.n = length(idx);

end

