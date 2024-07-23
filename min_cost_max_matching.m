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
    