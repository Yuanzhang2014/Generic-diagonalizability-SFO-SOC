% Check whether the matrix is generically diagonalizable. If yes,
% IfDiagonalizable=1; otherwise, IfDiagonalizable=0.
% Your Matlab version should support optimproblem and optimvar
%A = input('Please enter matrix A:\n');
function IfDiagonalizable=Check_Generic_Diagonalizability(a)
%For example: 
%a=zeros(5,5);
%a(1,2)=1;a(2,4)=1;a(4,2)=1;a(3,4)=1;a(5,4)=1;
%a(3,1)=1;
A=a;
IfDiagonalizable =0;
    A_sys = A';
    n=size(A,1);
    generic_rank = sprank(A_sys);
    
   % A1 = ones(n) - A_sys;
   % [min_matching1, min_totalcost1] = min_cost_max_matching(A1);
   % generic_rank = n - min_totalcost1;
    %generic_rank = hopcroftKarp(A_sys);
   % generic_rank = compute_generic_rank(A_sys);
%disp(['The generic rank of the structured matrix is: ' num2str(generic_rank)]); 

%璁＄畻MWMM(A)
    A_cost = zeros(n);
    for i = 1:n
        for j = 1:n
            if i == j && A_sys(i,j) == 0
                A_cost(i,j) = 1;
            end
            if A_sys(i,j) == 1
                A_cost(i,j) = 0;
            end
            if i~=j && A_sys(i,j) == 0
                A_cost(i,j) = 5000+1;
            end
        end
    end
    [min_matching, min_totalcost] = min_cost_max_matching(A_cost);
    %disp(['The MWMM of the structured matrix is: ' num2str(min_totalcost)]);
    if generic_rank == n - min_totalcost
    IfDiagonalizable =1;
     disp('Generically diagonalizable');  
    else
     disp('Not Generically diagonalizable');
    end
   
    
    %Solve the minimum weight maximum matching problem using linear programming
    % Your Matlab version should support optimvar.
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
end
    
