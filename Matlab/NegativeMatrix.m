function[answer] = NegativeMatrix(matrix)
    answer = true;
    [n,m] = size(matrix);
    if(n ~= m)
        error1 = "La matriz no es cuadrada"
        answer = False;
    end
    
    for i=1:n
        submatrix = matrix(1:i,1:i);
        module = mod(i,2)
        
        if(det(submatrix) == 0)
            error2 = 'La matriz no es negativa'
            answer = false;
            break;
        end
        
        if(module == 0 && det(submatrix) < 0)
            error2 = 'La matriz no es negativa'
            answer = false;
            break;
        elseif(module ~= 0 && det(submatrix) > 0 )
            error2 = 'La matriz no es negativa'
            answer = false;
            break;        
        end         
    end
end
