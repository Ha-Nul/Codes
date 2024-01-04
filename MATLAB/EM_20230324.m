x=linspace(-2,2,100);
f1=factorial(3)
f2=factorial(5)

y=sin(x);
y1=x;
y2=x-(1/3)*x.^3
y3=x-(1/f1)*x.^3+(1/f2)*x.^5

title('Sin(x)')
plot(x,y,'o')

hold on

plot(x,y1)
plot(x,y2,'*')
plot(x,y3,'m')

legend('sin(x)','y=x','y=x-x^3/3!','y=x-(x^3/3!)+(x^5/5!)')