function [nnodes,nelems,xNODE,xELEM] = X_1d_mesh(X)

nelems.x = length(X)-1; nnodes.x = length(X); % Number of elements & nodes

xELEM(:).nodes = zeros(1,nelems.x); % Allocate matricies
xELEM(:).jacobian = zeros(1,nelems.x);
for i = 1:nelems.x
    xELEM(i).nodes    = [i,i+1]; % nodes of each element
    xELEM(i).jacobian = 2/(X(i+1)-X(i)); % Jacobian of each element
end

% Set up node to element connectivity for periodic bcs
xNODE(:).elems = zeros(nelems.x, 2); % Allocate matrix
xNODE(1).elems = [ nelems.x,1 ];  % First node
for i = 2:nnodes.x-1
    xNODE(i).elems = [ i-1,i ];           % Internal nodes
end
xNODE(nnodes.x).elems = [ nelems.x,1 ];  % Last node

