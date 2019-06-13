function [answer] = SparseMatrix(A)
    [n,m] = size(A);
    contador = 0;
    for i=1:n
        for j=1:m
            if (A(i,j) == 0)
                contador = contador + 1;
            end
        end
    end
    
    if(contador > n)
        answer = true;
    else
        answer = false;
    end               
end

