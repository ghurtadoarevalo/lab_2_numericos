function [sol,error,costos]=cholesky(A, b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema

    [n,m]=size(A);
    costos = 1;
    for k = 1 : n
        A(k,k) = sqrt(A(k,k));
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        costos= costos+5;
        for j = k + 1 : n
            A(j:n,j) = A(j:n,j) - A(j,k)*A(j:n,k);
            costos = costos + 4;
        end
    end
    
    L = tril(A);
    z = inv(L)*b;
    sol = inv(L')*z;
    
    error = norm(eye(n)-inv(L*L')*A);
    costos = costos + 13;
end