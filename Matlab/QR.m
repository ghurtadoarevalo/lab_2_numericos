function [sol,error,costo] = QR(A, b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema
    costo = 0;
    [m, n] = size(A);
    R = zeros(n, n);
    V = A;
    Q=zeros(m, n);
    costo = costo + 4 ;
    for i =1:n
        R(i,i)= norm(V(:,i));
        Q(:,i)= V(:,i)/R(i,i);
        costo=costo+4;
        for j=i+1:m
           R(i,j)= (Q(:,i)')*V(:,j);
           V(:,j)=V(:,j) - R(i,j)*Q(:,i);
           costo = costo + 6 ;
        end
    end
    sol=inv(R)*Q'*b;
    error = norm(eye(m)-inv(Q*R)*A);
    costo= costo + 12;
end