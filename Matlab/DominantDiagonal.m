function[answer] = DominantDiagonal(matrix)
    answer = true;
    [n,m] = size(matrix);
    if(n ~= m)
        error1 = "La matriz no es cuadrada"
        answer = false;
    end
    
    for i=1:n
        diagonal = matrix(i,i);
        result = abs(sum(matrix(i, :)) - diagonal);
        if(diagonal <= result)
            answer = false;
        end
    end
end
