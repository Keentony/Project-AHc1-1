%Antonio Hernandez
%UH ID: 1168567
%Project AHc1-1
close, clc

%Stating constants 
A = 0; %lamda
PI = 4*atan(1);
D = dx^2/4;
N =50;

%Creating disretization of x and y with N+1 segments
%dx = dy
dx = (2*PI)/(N+1);


%Setting up corrdinate template for x and y
x = -PI:dx:PI;
y = -PI:dx:PI;

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
        F(ii,jj) = cos((x(ii)+2*PI)/2)*sin((y(ii)+PI)/2);
    end
end
%%
% Evaluating the Gauss-Seidel method
tol = 10^-6;
%Impementing max number of iterations that scales with changing N
it = 10*N^2;
error = 1+tol;
iterations = 0;
err = zeros(N+1,N);

while error > tol
    iterations = iterations+1;
    uold = U;
    for jj = 2:N+1
        for ii = 2:N+1
           %Updating within the boundaries except at Neumman BC
           U(ii,jj) = (0.25*(U(ii-1,jj)+U(ii+1,jj)+U(ii,jj-1)+U(ii,jj+1))-D*F(ii,jj));                      
        end
        %Applying Neummann BC
        U(N+2,jj) = (0.25*(2*U(N+1,jj)+U(N+2,jj-1)+U(N+2,jj+1))-D*F(N+2));
        for ii = 2:N+1
            %Calculating error from old u value to updated u value
            if U(ii,jj) ~= 0 
                err(ii,jj) = abs((U(ii,jj) - uold(ii,jj))/U(ii,jj)) * 100; 
            end 
        end
    end
   %Letting user know if you need more iterations to reach the tolerance  
   % of10^-6
   if iterations > it 
   disp('Max iteration was reached before tolerance was reached'); 
   break 
   end
   error = max(max(err));
   if error <= tol 
   fprintf('\n Tolerance reached at %d iterations \n',iterations);
   end
end
%%
%Plotting
[xx,yy]=meshgrid(x,y');
surf(yy,xx,U);
xlabel('x axis')
ylabel('y axis')
