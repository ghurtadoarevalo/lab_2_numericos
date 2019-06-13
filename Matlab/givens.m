function [sol,error,costo] = givens(A,b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema
    
    costo = 0;
    [n,m] = size(A);
     
    %if(n > 1000 && n <= 2000)
        %mitad de la matriz
    %    Nkeep = n/2;
    %    A=A(1:Nkeep,1:Nkeep);
    %    b=b(1:Nkeep);
    
    %elseif(n > 2000 && n <= 5000)
        %cuarto de la matriz
    %    Nkeep = n/4;
    %    A=A(1:Nkeep,1:Nkeep);
    %    b=b(1:Nkeep);
    %else
        %octavo de la matriz
    %    Nkeep = n/8;
    %    A=A(1:Nkeep,1:Nkeep);
    %    b=b(1:Nkeep);
    %end
    
    %[n,m] = size(A);
    
    Q = eye(n);
    R = A;
    costo = costo + 3;
    for i=1:m
        for k=i+1:n
            if (R(k,i) ~= 0)
                raiz = sqrt(R(k,i)^2 + R(i,i)^2);
                s = -R(k,i)/raiz;
                c = R(i,i)/raiz;
                G = eye(n); 
                
                G(k,k) = c;
                G(i,i) = c;
                
                G(i,k) = s;
                G(k,i) = -s;
                
                Q = Q*G; 
                R = G'*R; 
                costo = costo + 17;
            end
        end
    end

    Y = inv(Q)*b;  
    sol = inv(R)*Y;
    error = norm(eye(m)-inv(Q*R)*A);
    costo = costo + 13;
end
