%A naive iterative algorithm to determine a feasible solution to P1 for general systems
clear;
clc;

%Example 5 in the manuscript
%A = [0 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 1 0 1 0 0;0 1 0 0 0
%0];
%F = [0 1 1 1 0 0]

A = input('Please enter matrix A:\n');
F = input('Please enter matrix F:\n');    
n = length(A);
C = zeros(1,n);
h = length(F(:,1));
%Construct the set Xf
XF = zeros(1,n);
for i = 1 : h
    for j = 1 : n
        if F(i,j) ~= 0
            XF(j) = 1;
        end
    end
end
%Construct the output matrix
while gr_O(A,[C;F]) - gr_O(A,C) >=1
    C = [C;XF];
end
C = C(2:size(C,1),:); 

%Calculate the grank of the observability matrix. The randomized algorithm gives the right answer with probability one.
function [grank] = gr_O(A,C)
O = obsv(A,C);
p = size(O,1);
q = size(O,2);
grank = rank(O.*rand(p,q));
end
    