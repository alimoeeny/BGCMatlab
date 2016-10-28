function X = MatrixPermute(M, order)
%for square matrix M, re-oder rows and columnts according to order
%if order is shorter then length(M) returns a matix of length(order);
if isempty(order)
    X = [];
end
n = length(order);
for j = 1:n
    for k = 1:n;
        X(j,k) = M(order(j),order(k));
    end
end
