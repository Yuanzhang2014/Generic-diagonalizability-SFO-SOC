
function scc = stronglyConnectedComponents(adjMatrix)
    % Function to find strongly connected components in a directed graph
    % using Tarjan's algorithm
    %
    % Input:
    %   adjMatrix - NxN adjacency matrix of the directed graph
    %
    % Output:
    %   scc - cell array where each cell contains the nodes of one SCC

    n = size(adjMatrix, 1);
    index = 0;
    indices = -ones(n, 1);
    lowlink = -ones(n, 1);
    onStack = false(n, 1);
    stack = [];
    scc = {};
    
    for v = 1:n
        if indices(v) == -1
            [index, indices, lowlink, onStack, stack, scc] = strongconnect(v, index, indices, lowlink, onStack, stack, scc, adjMatrix);
        end
    end
end

function [index, indices, lowlink, onStack, stack, scc] = strongconnect(v, index, indices, lowlink, onStack, stack, scc, adjMatrix)
    % Helper function for Tarjan's algorithm
    index = index + 1;
    indices(v) = index;
    lowlink(v) = index;
    stack = [v stack];
    onStack(v) = true;
    
    for w = find(adjMatrix(v, :))
        if indices(w) == -1
            [index, indices, lowlink, onStack, stack, scc] = strongconnect(w, index, indices, lowlink, onStack, stack, scc, adjMatrix);
            lowlink(v) = min(lowlink(v), lowlink(w));
        elseif onStack(w)
            lowlink(v) = min(lowlink(v), indices(w));
        end
    end
    
    if lowlink(v) == indices(v)
        new_scc = [];
        while true
            w = stack(1);
            stack(1) = [];
            onStack(w) = false;
            new_scc = [new_scc w];
            if w == v
                break;
            end
        end
        scc{end+1} = new_scc;
    end
end
