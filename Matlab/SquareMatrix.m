function [answer] = SquareMatrix(A)
    [n,m] = size(A);
    if( n == m)
        answer = true;
    else
        answer = false;
    end
end