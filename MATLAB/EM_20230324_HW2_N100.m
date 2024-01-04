clearvars;

a = 20;
b = 10;
Na = 100;
Nb = 100;
V01 = 1.0;
V02 = 2.0;
V03 = 3.0;
V04 = 1.0;
[X,Y] = meshgrid(0:a/Na:a,0:b/Nb:b);

Nk = 10;
for k = 1:Nk
    n = 2*k-1;
    V1(:,:,k) = 4*V01/pi*sin(n*pi*X/a).*sinh(n*pi*Y/a)./(n*sinh(n*pi*b/a));
end

for k = 1:Nk
    n = 2*k-1;
    V2(:,:,k) = 4*V02/pi*sin(n*pi*X/a).*sinh(n*pi*Y/a)./(n*sinh(n*pi*b/a));
end

for k = 1:Nk
    n = 2*k-1;
    V3(:,:,k) = 4*V03/pi*sin(n*pi*X/a).*sinh(n*pi*Y/a)./(n*sinh(n*pi*b/a));
end

for k = 1:Nk
    n = 2*k-1;
    V4(:,:,k) = 4*V04/pi*sin(n*pi*X/a).*sinh(n*pi*Y/a)./(n*sinh(n*pi*b/a));
end

V1_Tot=sum(V1,3);
V2_Tot=sum(V2,3);
V3_Tot=sum(V3,3);
V4_Tot=sum(V4,3);

V_total=V1_Tot+V2_Tot+V3_Tot+V4_Tot

surf(X,Y,V_total)
xlim=[0,20];
ylim=[0,10];
zlim=[0,7];

% In the case of N=100, compare with the case of N=10, the total summation value approximates to
% Z=(const). We can assume if value of N goes larger, then we can get more
% correct value for Z components.