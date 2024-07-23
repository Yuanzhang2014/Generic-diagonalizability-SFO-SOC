% Algorithm 1 : A weighted maximum matching based algorithm for P1 in generically diagonalizable systems
clear;
clc;
A = input('Please enter matrix A:\n');
F = input('Please enter matrix F:\n');    
n = length(A(1,:));
h = length(F(:,1));
%Construct the set Xf
XF = zeros(n,1);
flag=Check_Generic_Diagonalizability(A);
if flag==0
    disp('The input system is not generically diagonalizable. The algorithm can only give a lower bound.')
end
for i = 1 : h
    for j = 1 : n
        if F(i,j) ~= 0
            XF(j) = 1;
        end
    end
end
%Construct the weight matrix
A_cost = zeros(n);
for i = 1 : n
    for j = 1 : n
        if A(i,j) ~= 0 && XF(j) == 1
            A_cost(i,j) = 1;
        elseif A(i,j) ~= 0 && XF(j) == 0
            A_cost(i,j) = 0;
        else A_cost(i,j) = 2 * n+2;
        end
    end
end
A_cost = A_cost';
%Find the solution for the minimum weight maximum matching
[min_matching] = min_cost_max_matching(A_cost);
for i = 1 : n
    for j = 1 : n
        if min_matching(i,j) == 1 && A(j,i) == 0
            min_matching(i,j) = 0;
        end
    end
end
%Construct the set Xf-unmatched
X_m = zeros(n,1);
XF_u = zeros(n,1);
for i = 1 : n
    for j = 1 : n
        if min_matching(i,j) == 1
            X_m(i) = 1;
        end
    end
end
for i = 1 : n
    if X_m(i) == 0 && XF(i) == 1
        XF_u(i) = 1;
    end
end
%Construct the output matrix
C = zeros(sum(XF_u),n);
flag = 1;
for i = 1 : n
    if XF_u(i) == 1
        C(flag,i) = 1;
        flag = flag + 1;
    end
end
for i = 1 : n
    if X_m(i) == 1 && XF(i) == 1
        C(1,i) = 1;
    end
end

%Solve the minimum weight maximum matching problem using linear programming
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
    
