function[answer] = SemiPositiveMatrix(matrix)
    answer = true;
    [n,m] = size(matrix);
    if(n ~= m)
        error1 = "La matriz no es cuadrada"
        answer = False;
    end
    
    for i=1:n
        submatrix = matrix(1:i,1:i);
        if(det(submatrix) <= 0)
            error2 = 'La matriz no es semi positiva'
            answer = false;
            break;
        end
    end
end
