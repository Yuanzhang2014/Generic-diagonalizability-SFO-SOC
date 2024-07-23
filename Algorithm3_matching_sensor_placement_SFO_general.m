%Algorithm 3 : A weighted maximum matching based algorithm for a feasible solution to P1 for general systems
clear;
clc;

%Example 5 in the manuscript
%A = [0 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 1 0 1 0 0;0 1 0 0 0
%0];
%F = [0 1 1 1 0 0]

A = input('Please enter matrix A:\n');
F = input('Please enter matrix F:\n');    
n = length(A(1,:));
h = length(F(:,1));
XF = zeros(n,1);
%Construct the set Xf
for i = 1 : h
    for j = 1 : n
        if F(i,j) ~= 0
            XF(j) = 1;
        end
    end
end
F_set = null(1);
for i = 1 : n
    if XF(i) == 1
        F_set = [F_set, i];
    end
end
%Construct the output reachable set WF
Y_reachable = findNodesReaching(A', F_set);
%Construct the weight matrix
q = sum(XF);
%The weight of Exx
A_cost1 = zeros(n);
for i = 1 : n
    for j = 1 : n
        if A(i,j) ~= 0
            if find(Y_reachable == i)
                A_cost1(i,j) = -(q + 1);
            else
                A_cost1(i,j) = 0;
            end
        elseif i == j
            A_cost1(i,j) = 0;
        else
            A_cost1(i,j) = 2 * n * n + 2 * n + 1;%Assign a large weight to represent no edge
        end
    end
end
%The weight of Exy
A_cost2 = (2 * n * n + 2 * n + 1) .* ones(q,n);
flag1 = 1;
for i = 1 : n
    if find(F_set == i)
        A_cost2(flag1,i) = -q;
        flag1 = flag1 + 1;
    end
end
%The weight of Eyx
A_cost3 = zeros(n,q);
%The weight of Eyy
A_cost4 = (2 * n * n + 2 * n + 1) .* ones(q) - (2 * n * n + 2 * n + 1) .* eye(q);
A_cost = [A_cost1 A_cost3; A_cost2, A_cost4];
%Find the solution for the minimum weight maximum matching
A_cost = A_cost';
[min_matching] = min_cost_max_matching(A_cost);
%Construct the set Xh
XH = zeros(1, n);
for i = n + 1 : n + q
    for j = 1 : n
        if min_matching(j,i) == 1
            XH(j) = 1;
        end
    end
end
%Construct the output matrix
C = zeros(sum(XH),n);
flag = 1;
for i = 1 : n
    if XH(i) == 1
        C(flag,i) = 1;
        flag = flag + 1;
    end
end
%Construct the H-reachable set H_reachable
H_set = null(1);
for i = 1 : n
    if XH(i) == 1
        H_set = [H_set, i];
    end
end
H_reachable = findNodesReaching(A', H_set);
for i = 1 : n
    if XH(i) ~= 1 && XF(i) == 1 && ~find (H_reachable == i)
        C(1,i) = 1;
    end
end

function [matching, total_cost] = min_cost_max_matching(cost_matrix)
    n = size(cost_matrix, 1);

    % Define variables
    x = optimvar('x', n, n, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);

    % Construct the objective function
    obj = sum(sum(cost_matrix .* x));

    % Construct the constraints
    constr = [sum(x, 1) == 1;
              sum(x, 2)' == 1];

    % Create the problem object
    prob = optimproblem('Objective', obj, 'Constraints', constr);

    % Solve the problem
    [sol, fval, exitflag] = solve(prob);

    % Analyze the results.
    if exitflag == 1
        matching = sol.x;
        total_cost = fval;
    else
        error('Solution failed');
    end
end
    
function reachableNodes = findNodesReaching(adjMatrix, targetNodes)
    % Obtain the number of nodes
    numNodes = size(adjMatrix, 1);
    
    % Inverse adjacency matrix
    reverseAdjMatrix = adjMatrix';
    
    % Record the nodes that can reach any target node
    visited = false(1, numNodes);
    
    % Use a stack to perform DFS, starting from each target node.
    for targetNode = targetNodes
        stack = [targetNode]; 
        while ~isempty(stack)
            currentNode = stack(end);
            stack(end) = [];
            if ~visited(currentNode)
                visited(currentNode) = true;
                neighbors = find(reverseAdjMatrix(currentNode, :));
                for neighbor = neighbors
                    if ~visited(neighbor)
                        stack(end + 1) = neighbor;
                    end
                end
            end
        end
    end
    
    % Return the nodes that can reach any target node
    reachableNodes = find(visited);
end