%Antonio Hernandez
%UH ID: 1168567
%Project AHc1-1
close, clc
%Stating constants 
A = 1; %lamda
PI = 4*atan(1);
ax = -PI;
ay = -PI; 
bx = PI;
by = PI;
T = 10; %Max time value

%Creating disretization of x and y with N+1 segments
N = 50;
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

%Setting up the three Direlecht boundary conditions 
%Stating q and w outside loop to increase optimization
U = zeros(N+2,N+2);
q = -4*(PI)^3;
w = 4*(PI)^2*cos(-PI)+4*(PI)^3;
for ii = 1:N+2
    U(ii,1)=x(ii)*(PI-x(ii))^2;
    U(ii,N+2) = (PI-x(ii))^2*cos(x(ii));
    U(1,ii) = q+((x(ii)+PI)/(2*PI))*w;
end

%Setting the values of F(x,y) from -PI to PI for x and y
F = zeros(N+2,N+2);
for ii = 1:N+2
    for jj = 1:N+2
        F(ii,jj) = cos((PI/2)*(2*((x(ii)-ax)/(bx-ax))+1))*sin(PI*((y(jj)-ay)/(by-ay)));
    end
end

%Extracting the desired values of F(x,y) where x is up to N+2 because
%of the Neumann condition and y is up to N+1
F = F(2:N+2,2:N+1);
F = dx^2*F;

%Adding the three Direlecht boundary conditions to the F(x,y) matrix
F(:,1)=F(:,1)-U(2:N+2,1);
F(:,N)=F(:,N)-U(2:N+2,N);
F(1,:)=F(1,:)-U(1,2:N+1);

%Reshaping the Fij and uij matrices into their proper orientation
F = reshape(F,(N+1)*N,1);
u = U(2:N+2,2:N+1);
u = reshape(u,(N+1)*N,1);

%Setting the diagonal matrix for the constants of U with the Neumann
%condition where the constant is 2 at every N+1 row and Nth collumn of H
%Note: H and I are the building blocks that lead to the overall matrix J
H = diag(ones(1,N),-1)+diag(D*ones(1,N+1))+diag(ones(1,N),1);
H(N+1,N) = 2;
I = diag(ones(1,N+1));

%Bringing the H and I matrices together in a 'checkerd' formation to create
%the overall matrix J which includes Neumann condition
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
