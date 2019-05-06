function [nr_edge_adj, mat_out] = determineAdjacentEdges(edge_set, vrtx )

nr_edge_adj=0;

%% Description 
% Find edges in the set edge_set that are adjacent to the vertex vrtx
% Input.
%   edge_set: the set of edges
%   vrtx: the reference vertex
% Output.
%   nr_edge_adj: the number of adjacent edges
%   mat_out: matrix of adjacent edges

[nr_edges,~]=size(edge_set);
mat_out=zeros(1,nr_edges);
for ii=1:nr_edges
   if edge_set(ii,1)==vrtx || edge_set(ii,2)==vrtx
      nr_edge_adj=nr_edge_adj+1;
      mat_out(1,nr_edge_adj)=ii;    
   end
    
end
mat_out=setdiff(mat_out, 0);

end

