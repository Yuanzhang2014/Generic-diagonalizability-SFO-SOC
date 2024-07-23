clear;
clc;

q = 10;
Q = zeros (2,10);
Q(1,1) = 100;
Q(1,2) = 150;
Q(1,3) = 200;
Q(1,4) = 250;
Q(1,5) = 300;
Q(1,6) = 800;
Q(1,7) = 850;
Q(1,8) = 900;
Q(1,9) = 950;
Q(1,10) = 1000;
%Q(2,1) = 0.025;
%Q(2,2) = 0.05;
%Q(2,3) = 0.075;
%Q(2,4) = 0.1;
K = zeros(2,length(Q(1,:)));

 for nl = 1:length(Q(1,:))
    flag = 0;
    s = log(Q(1,nl))/Q(1,nl);
    for i = 1:q
    A = ER_network(Q(1,nl),s); % 生成 ER 随机网络
% 可视化邻接矩�?
%G = digraph(A);
%plot(G);

% 计算泛秩
    A_sys = A';
    generic_rank = rank(A_sys .* rand(Q(1,nl)));
   % A1 = ones(Q(1,nl)) - A_sys;
   % [min_matching1, min_totalcost1] = min_cost_max_matching(A1);
   % generic_rank = Q(1,nl) - min_totalcost1;
    %generic_rank = hopcroftKarp(A_sys);
   % generic_rank = compute_generic_rank(A_sys);
%disp(['The generic rank of the structured matrix is: ' num2str(generic_rank)]); 

%计算MWMM(A)
    A_cost = zeros(Q(1,nl));
    for i = 1:Q(1,nl)
        for j = 1:Q(1,nl)
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
    if generic_rank == Q(1,nl) - min_totalcost
    flag = flag +1;
    end
    end
     K(2,nl) = flag/q;
 end
K(1,:) = Q(1,:);
save('C:\Users\HUAWEI\Desktop/data.mat', 'K'); 
%plot(Q(1,:),K(1,:),'*-','Color','g');
%hold on;
%plot(Q(1,:),K(2,:),'o-','Color','r');
%plot(Q(1,:),K(3,:),'+-','Color','b');
%plot(Q(1,:),K(4,:),'x-','Color','k');
plot(K(1,:),K(2,:),'s-','Color','k');
%legend('p=0.025','p=0.05','p=0.075','p=0.1','p=log(n)/n');
legend('p=log(n)/n');
title('ER networks q=100');
xlabel('n');
ylabel('Proportion')

    
