function x = graph1(fil,col)
%matriz random 10x10
A = rand(10,10);
%matriz random 1x10
b = rand(1,10);
%matriz transpuesta
b_transpuesto = b';
x = A*b_transpuesto

%grafico de c
%plot(c)

figure1 = figure;
createfigure2(A,figure1);
end