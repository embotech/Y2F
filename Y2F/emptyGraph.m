function G = emptyGraph()
%EMPTYGRAPH returns an empty graph G
% Graphs are stored as structs with the following fields:
%   .vertices   cell array of vertex labels (lists of variable indices)
%   .adjMatrix  adjacency matrix
%   .n          number of vertices

G.vertices = {};
G.adjMatrix = [];
G.n = 0;

end

