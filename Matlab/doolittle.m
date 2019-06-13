function [sol,error,costos] = doolittle(A,b)

    % sol es vector solución
    % error es el error final
    % A matriz del sistema
    % b vector del sistema
    costos = 0;
    [n, m]=size(A);
    L=zeros(n);
    U=zeros(n);
    L(1,1)=1;
    U(1,1)=A(1,1);
    costos = costos + 5;
    for i=2:n
        L(i,i)=1;
        U(1,i)=A(1,i);
        L(i,1)=A(i,1)/U(1,1);
        costos = costos + 3;
    end
    for j=2:n
        for i=j:n
            sumal=0;
            sumau=0;
            costos = costos + 2;
            for k=1:j-1
                if (U(k,i)~=0) && (L(j,k)~=0)
                    sumal=sumal+U(k,i)*L(j,k);
                    costos = costos + 7;
                end
                if (U(k,j)~=0)&&(L(i,k)~=0)&&(i~=j)
                   sumau=sumau+U(k,j)*L(i,k);
                   costos = costos + 8;
                end
            end
                U(j,i)=A(j,i)-sumal;
                costos = costos + 2;
            if (j<n)&&(i>j)
                L(i,j)=(A(i,j)-sumau)/U(j,j);
                costos = costos + 6;
            end
        end
    end
        
    z = inv(L)*b;
    sol = inv(U)*z;
    error = norm(eye(n)-inv(L*U)*A);  
    costos = costos + 12;
end
