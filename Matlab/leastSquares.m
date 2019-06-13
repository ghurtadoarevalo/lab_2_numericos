function [error,sol] = leastSquares(A, b)
    [n,m] = size(A);
    sol = inv(A'*A)*A'*b;
    error = norm(eye(n,m)-(A'*A)*A');
end