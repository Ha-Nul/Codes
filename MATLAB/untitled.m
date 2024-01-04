clearvars;

a = 20;
b = 10;
Na = 100;
Nb = 100;
V0 = 1.0;
[X,Y] = meshgrid(0:a/Na:a,0:b/Nb:b);

Nk = 100;
for k = 1:Nk
    n = 2*k-1;
    V(:,:,k) = 4*V0/pi*sin(n*pi*X/a).*sinh(n*pi*Y/a)./(n*sinh(n*pi*b/a));
end

V_total = sum(V,3);
    
surf(X,Y,V_total)