function [maxFlow, minCost, pathMatrices] = minCostMaxFlow(capacity, cost, source, sink)
    % Initialize parameters
    n = size(capacity, 1);  % Number of nodes
    flow = zeros(n);  % Flow matrix
    maxFlow = 0;  % Total maximum flow
    minCost = 0;  % Total minimum cost
    pathMatrices = {};  % Save all augmented path matrices
    % Loop until there are no augmented paths
    while true
        % Step 1: Use the Bellman-Ford algorithm to find the shortest path from the source to the sink
        [dist, parent] = bellmanFord(n, capacity, cost, flow, source);
        
        % If there is no path from the source to the sink, exit the loop
        if dist(sink) == Inf
            break;
        end
        
        % 步骤2：Step 2: Determine the maximum flow of the found path
        increment = Inf;
        v = sink;
        while v ~= source
            u = parent(v);
            increment = min(increment, capacity(u, v) - flow(u, v));
            v = u;
        end
        
        % Record the current path's
        pathMatrix = zeros(n);
        v = sink;
        while v ~= source
            u = parent(v);
            pathMatrix(u, v) = 1;  % Mark the edges in the path
            v = u;
        end
        pathMatrices{end + 1} = pathMatrix;  % Save the path matrix
        
        % 步骤3：Update flow and cost
        v = sink;
        while v ~= source
            u = parent(v);
            flow(u, v) = flow(u, v) + increment;
            flow(v, u) = flow(v, u) - increment;
            minCost = minCost + increment * cost(u, v);
            v = u;
        end
        
        maxFlow = maxFlow + increment;
    end
end

function [dist, parent] = bellmanFord(n, capacity, cost, flow, source)
    dist = Inf * ones(1, n);  % Distance from the source to each node
    dist(source) = 0;
    parent = -1 * ones(1, n);  % Parent node of each node in the path

    % Relax at most n-1 times
    for i = 1:n-1
        for u = 1:n
            for v = 1:n
                if capacity(u, v) > flow(u, v) && dist(v) > dist(u) + cost(u, v)
                    dist(v) = dist(u) + cost(u, v);
                    parent(v) = u;
                end
            end
        end
    end
end
