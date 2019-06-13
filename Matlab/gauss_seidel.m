function [sol,errors, costos] = gauss_seidel(A, b)
    % sol es vector solución
    % err es el error final
    % it es el numero de iteraciones final
    % A matriz del sistema
    % b vector del sistema
    % maxiter es numero máximo de iteraciones
    % epsilon es la cota del error
    costos = 0;
    n=length(b);
    sol=zeros(1,n)';

    error = norm(A*sol - b);

    tolerance = 7.326146437847726e-16;
    iterations = 0;

    if(n <= 1000)
        maxiter = 100;
    elseif(n > 1000 && n<= 4000)
        maxiter = 500;
    else
        maxiter = 1000;
    end

    
    errors = [];
    errors = [errors, error];
    costos = costos + 6;
    while error > tolerance && iterations < maxiter
        sol_old = sol;  
        costos = costos + 1;
        for i=1:n
            S=0;
            costos = costos + 1;
            for j=1:i-1
                S=S+A(i,j)*sol(j);
                costos = costos + 1;
            end

            for j=i+1:n
                S=S+A(i,j)*sol_old(j);
                costos = costos + 1;
            end

            sol(i)=(1/A(i,i))*(b(i)-S);
        end

        iterations=iterations+1;
        error=norm(sol_old-sol);
        errors = [errors, error];
        costos = costos + 3;
    end
end