%newton rapson necesita una solucion, tolerancia, cantidad de iteraciones,
%la funcion y su primera derivada
function [xVector,errores,iteracionesCount] = newtonRapson(x0,tolerancia, iteraciones, fx, dfx)
x = x0;
errores = [];
xVector = [];
for i=1:iteraciones
    iteracionesCount = i;
    error = abs(fx(x))
    if abs(fx(x)) < tolerancia || abs(fx(x)) == 0  
        break;
    end
    if errores < tolerancia
        break;
    end
    xVector = [xVector, x]
    errores = [errores, error]
    
    x = x - fx(x)/dfx(x);

end
end

%matriz random 10x10
%A = rand(10,10);
%matriz random 1x10
%b = rand(1,10);
%matriz transpuesta
%b_transpuesto = b';
%x = A*b_transpuesto

%grafico de c
%plot(c)

%figure1 = figure;
%createfigure2(A,figure1);