function[answer] = SimetricMatrix(matrix)
    answer = true;
    [n,m] = size(matrix);
    if(n ~= m)
        answer = False;
    end
    
    for i=1:n
        for j=1:m
            if(matrix(i,j) ~= matrix(j,i))
                answer = false;
                break;
            end
        end
    end
end
