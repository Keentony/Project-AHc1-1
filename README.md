%Antonio Hernandez
%UH ID: 1168567
%Project AHc1-1

%Stating constants 
A = 1; %lamda
PI = 4*atan(1);
ax = -PI;
ay = -PI; 
bx = PI;
by = PI;

%Creating disretization of x and y with N+1 segments
N = 10;
dx = (bx-ax)/(N+1);
dy = (by-ay)/(N+1);

%Constant for main diagonal of U(i,j)
D = -4+(dx^2*A);

%Setting up corrdinate template for x and y
x = zeros(N+2,1);
x(1)=-PI;
for ii = 1:N+1
    x(ii+1) = -PI+(ii*dx);
end
y = zeros(N+2,1);
y(1)=-PI;
for jj = 1:N+1
    y(jj+1) = -PI+(jj*dy);
end

%Setting the (x,y)
F = zeros(N+2,N+2);
for ii = 1:N+2
    for jj = 1:N+2
        F(ii,jj) = cos((PI/2)*(2*((ii-ax)/(bx-ax))+1))*sin(PI*((jj-ay)/(by-ay)));
    end
end
