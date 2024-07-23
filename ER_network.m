function A = ER_network(n,p)
% ER_network - Generate an ER random network
% n - Number of network nodes
% p - Probability of edges
 
A = zeros(n); % Initialize the adjacency matrix as an all-zero matrix
 
for i = 1:n
    for j = 1:n
        if rand() <= p %If the random number is less than or equal to the probability p, add an edge between nodes i and j
            A(i,j) = 1;
        end
    end
end

end 
