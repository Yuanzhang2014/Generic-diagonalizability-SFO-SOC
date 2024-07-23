function A = ER_network(n,p)
% ER_network - 生成 ER 随机网络
% n - 网络节点数
% p - 边的概率
 
A = zeros(n); % 初始化邻接矩阵为全零矩阵
 
for i = 1:n
    for j = 1:n
        if rand() <= p % 如果随机数小于等于概率 p，则在节点 i 和 j 之间添加一条边
            A(i,j) = 1;
        end
    end
end

end %ER网络生成