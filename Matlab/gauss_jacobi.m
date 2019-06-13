function [sol,errors, costos] = gauss_jacobi(A, b)  
    %Si A es una matriz de diagonal estrictamente dominante por filas,
    %El método converge.
    costos = 0;
  
    epsilon = 7.326146437847726e-16;
    
    % sol es vector solución
    % error es el error final
    % it es el numero de iteraciones final
    % A matriz del sistema
    % b vector del sistema
    % maxiter es numero máximo de iteraciones
    % epsilon es la cota del error
    
    n=length(b);
    
    if(n <= 1000)
        maxiter = 100;
    elseif(n > 1000 && n<= 4000)
        maxiter = 500;
    else
        maxiter = 1000;
    end

        
    
    iterations=0;
    sol_x=zeros(1,n)';
    sol=zeros(1,n);

    error = norm(A*sol_x - b);
    errors = [];
    errors = [errors, error];
    costos = costos + 9;
    while error>epsilon && iterations < maxiter
        for i=1:n
           S=0;
           costos = costos + 1;
           for j=1:n
               if i~=j
                S=S+A(i,j)*sol_x(j);
                costos = costos + 1;
               end
           end
           sol(i)=(b(i)-S)/A(i,i);
           costos = costos + 1;
        end
        iterations=iterations+1;
        error=norm(sol_x-sol);
        errors = [errors, error];
        sol_x=sol;
        costos = costos + 4;
    end
    sol = sol';
    costos = costos + 1;
end