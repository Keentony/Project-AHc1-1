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
N = 3;
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

%Setting the values of the function from -PI to PI for x and y
F = zeros(N+2,N+2);
for ii = 1:N+2
    for jj = 1:N+2
        F(ii,jj) = cos((PI/2)*(2*((x(ii)-ax)/(bx-ax))+1))*sin(PI*((y(jj)-ay)/(by-ay)));
    end
end

%Extracting the desired values of the function where x is up to N+2 because
%of the Neumann condition 
F = F(2:N+2,2:N+1);

%Setting the diagonal matrix for the constants of U with the Neumann
%condition where the constant is 2 at every N+1 row and Nth collumn of H
%Note: H and I are the building blocks that lead to the overall matrix
H = diag(ones(1,N),-1)+diag(D*ones(1,N+1))+diag(ones(1,N),1);
H(N+1,N) = 2;
I = diag(ones(1,N+1));
%Bringing the H and I matrices together in a 'checkerd' formation to create
%the overall matrix
J = sparse((N+1)*N,(N+1)*N);
for ii = 1:N
    J(1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1),1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1))=H;
end

for jj = 1:N-1
    J(jj*(N+1)+1:jj*(N+1)+(N+1),1+(jj-1)*(N+1):(N+1)+(jj-1)*(N+1))=I;
end

for kk = 1:N-1
    J(1+(kk-1)*(N+1):(N+1)+(kk-1)*(N+1),kk*(N+1)+1:kk*(N+1)+(N+1))=I;
end
