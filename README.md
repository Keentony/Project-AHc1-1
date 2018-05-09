%Antonio Hernandez
%UH ID: 1168567
%Project AHc1-1
close, clc


%Stating constants 
% With lambda = 0 we are solving Poisson's 2D equation
PI = 4*atan(1);
% ax, bx, ay, and by were left as is for the sake of modularity
ax = -PI;
bx = PI;
ay = -PI;
by = PI;
%Number of Segments
N = 75;

%Creating discretization of x and y with N+1 segments
%In this case solving for Poisson with dx = dy
dx = (bx-ax)/(N+1);
D = dx^2/4;

%Setting up corrdinate template for x and y
x = ax:dx:bx;
y = ay:dx:by;

%Setting up the three Direlecht boundary conditions 
%Stating q and w outside the loop to increase optimization
U = zeros(N+2,N+2);
gb = ax*(bx-ax)^2;
fb_gb = ((bx-ax)^2*cos(PI*ax/bx))-(ax*(bx-ax)^2);
for ii = 1:N+2
    U(ii,1)=x(ii)*(bx-x(ii))^2;
    U(ii,N+2) = (bx-x(ii))^2*cos((PI*x(ii))/bx);
    U(1,ii) = gb+(((y(ii)-ay))/(by-ay))*fb_gb;
end
%Setting the values of F(x,y) from -PI to PI for x and y
F = zeros(N+2,N+2);
PI_2=PI/2;
for ii = 1:N+2
    for jj = 1:N+2
        F(ii,jj) = cos(PI*(x(ii)-ax)/(bx-ax)+PI_2)*sin(PI*((y(jj)-ay)/by-ay));
    end
end

%%
% Evaluating the Gauss-Seidel and SOR method
% Gauss-Seidel is when w=1 and SOR can have a range of 1<w<2
tol = 10^-6;
%Impementing max number of iterations that scales with changing N 
%Higher N requires higer amounts of iterations to reach the tolerance
it = 9*N^2;
error = 1+tol;
iterations = 0;
err = zeros(N+1,N);
tic
while error > tol
    iterations = iterations+1;
    uold = U;
    w = 1.7; %For SOR w = 1.7 is nominal for different N values
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
                %Applying SOR 
                U(ii,jj) = w*U(ii,jj)+(1-w)*uold(ii,jj);
                err(ii,jj) = abs((U(ii,jj) - uold(ii,jj))/U(ii,jj)) * 100; 
            end 
        end
    end
   %Letting user know if you need more iterations to reach the tolerance  
   % of 10^-6
   if iterations > it 
   disp('Max iteration was reached before tolerance was reached'); 
   break 
   end
   error = max(max(err));
   if error <= tol 
   fprintf('Tolerance reached at %d iterations \n',iterations);
   end
end
toc

%%
%Plotting
[xx,yy]=meshgrid(x,y');
surf(yy,xx,U);
xlabel('x axis')
ylabel('y axis')
