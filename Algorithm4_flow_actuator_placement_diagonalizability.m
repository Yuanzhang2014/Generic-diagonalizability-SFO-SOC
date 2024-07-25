clear;
clc;
%Example
% A = [0 1 0 0 0;
%     0 0 0 1 0;
%     0 0 0 1 0;
%     0 1 0 0 0;
%     0 0 0 1 0];
%  C = [1 0 1 0 0;
%      0 1 0 0 0;
%      0 0 0 0 1];
A = input('Please enter matrix A:\n');
C = input('Please enter matrix C:\n');
flag=Check_Generic_Diagonalizability(A);
if flag==0
    disp('The input system is not generically diagonalizable. The algorithm can only give a lower bound.')
else
    disp('The input system is generically diagonalizable.')
end
n = size(A,1);
q = size(C,1);
for j = 1 : n
    for i = 1 : n
        if A(i,j) ~= 0
            A(i,j) = 1;
        end
    end
    for i = 1 : q
        if C(i,j) ~= 0
            C(i,j) = 1;
        end
    end
end
%Construct the flow network adjacency matrix
G = zeros(6 * n + 2 * q + 2);
G(n + 1 : 2 * n, 1 : n) = eye(n);
G(2 * n + 1 : 3 * n, n + 1 : 2 * n) = A;
G(3 * n + 1 : 4 * n, 2 * n + 1 : 3 * n) = eye(n);
G(6 * n + 1 : 6 * n + q, 3 * n + 1 : 4 * n) = C;
G(5 * n + 1 : 6 * n, 4 * n + 1 : 5 * n) = eye(n);
G(2 * n + 1 : 3 * n, 5 * n + 1 : 6 * n) = eye(n);
G(6 * n + q + 1 : 6 * n + 2 * q, 6 * n + 1: 6 * n + q) = eye(q);
G(1 : n, 6 * n + 2 * q + 1) = ones(n, 1);
G(4 * n + 1 : 5 * n, 6 * n + 2 * q + 1) = ones(n, 1);
G(6 * n + 2 * q + 2, 6 * n + q + 1 : 6 * n + 2 * q) = ones(1, q);
%Construct the cost matrix
Cost_matrix = zeros(6 * n + 2 * q + 2);
Cost_matrix(2 * n + 1 : 3 * n, 5 * n + 1 : 6 * n) = eye(n);
%Determine the source and sink nodes
source = 6 * n + 2 * q + 1;
sink = 6 * n + 2 * q + 2;
%Solve the maximum flow problem
[maxFlow, minCost, pathMatrices] = minCostMaxFlow(G', Cost_matrix', source, sink);
TotalFlow = zeros(6 * n + 2 * q + 2);
for i = 1 : length(pathMatrices)
    TotalFlow = TotalFlow + pathMatrices{i};
end
%Construct Xf1
X = zeros(1,n);
flag = 0;
Xf1 = null(1);
Xf2 = null(1);
for i = 5 * n + 1 : 6 * n
    for j = 2 * n + 1 : 3 * n
        if TotalFlow(i,j) ~= 0
            X(j - 2 * n) = 2;
            flag = flag + 1;
            Xf1 = [Xf1, j - 2 * n];
        end
    end
end
%Construct Xf2
for i = n + 1 : 2 * n
    for j = 2 * n + 1 : 3 * n
        if TotalFlow(i,j) ~= 0
            X(j - 2 * n) = 1;
            Xf2 = [Xf2, j - 2 * n];
        end
    end
end
B = zeros(flag, n);
%Strongly connected decomposition
scc = stronglyConnectedComponents(A');
%Construct B
for i =1 : n
    if X(i) == 2
        B(flag, i) = 1;
        flag = flag - 1;
    end
end
for i = 1 : length(scc)
    if ~isempty(intersect(scc{i},Xf2)) && isempty(intersect(scc{i},Xf1))
        B(1,scc{i}(length(scc{i}))) = 1;
    end
end
% Output the solution. Notice that the solution is not unique. 
disp('The solution is:')
disp(B)


