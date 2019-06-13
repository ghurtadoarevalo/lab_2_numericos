% Metodo de sustitucion regresiva para
% sistemas de ecuaciones lineales
function [x] = sustreg(A, b)
    % ¿Es válido el sistema?
    [M,N] = size(A);
    tam = length(b);
    if (M ~= N)||(tam ~= M)
        error('Las dimensiones del sistema no son válidas');
    end
    % Vector de soluciones
    x = zeros(tam,1);
    x(tam) = b(tam)/A(tam,tam);
    for k = tam-1:-1:1
        suma = x(k+1:tam)'*A(k,k+1:tam)';
        x(k) = (b(k)-suma)/A(k,k);
    end
end
