function w=poissondoubleP(xl,xr,yb,yt,M,N,L_start,L_stop,R_start,R_stop,p,p2, K, H, plot)
% Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps
% Output: matrix w holding solution values
% Example use: w = poisson(0,6,0,6,50,50,2,4,1,5,2,7,1.68,0.005,true)
%Power Variables
delta = 0.1;
L = L_stop - L_start;
L2 = R_stop-R_start;
f=@(x,y) 0; % define input function data
P=@(y) -(p/(L*delta*K));
P2 = @(y) -(p2/(L2*delta*K));
g1=@(x) 0; % define boundary values
g2=@(x) 0; % Example 8.8 is shown
g3=@(y) 0;
g4=@(y) 0;

%Constants
M = M - 1;
N = N - 1;
m=M+1;n=N+1; mn=m*n;
h=(xr-xl)/M;
h2=h^2;
k=(yt-yb)/N;
k2=k^2;
x=xl+(0:M)*h; % set mesh values
y=yb+(0:N)*k;

A=zeros(mn,mn);b=zeros(mn,1); %Initialize A and b

for i=2:m-1 % interior points
    for j=2:n-1
        A(i+(j-1)*m,i+(j-1)*m)=-2*((1/h2)+(1/k2)+(H/(K*delta)));
        A(i+(j-1)*m,i+1+(j-1)*m)=1/h2;
        A(i+(j-1)*m,(i-1)+(j-1)*m)=1/h2;
        A(i+(j-1)*m,i+(j-2)*m)=1/k2;
        A(i+(j-1)*m,i+j*m)=1/k2;
        b(i+(j-1)*m)=f(x(i),y(j));
    end
end
for i=2:m-1 % bottom and top boundary points
    j=1;
    A(i+(j-1)*m,i+(j-1)*m)= ((-3/(2*k))+(H/K));
    A(i+(j-1)*m,i+j*m) = (2/k);
    A(i+(j-1)*m,i+(j+1)*m) = -(1/(2*k));
    b(i+(j-1)*m)=g1(x(i));

    j=n;
    A(i+(j-1)*m,i+(j-1)*m)= ((-3/(2*k))+(H/K));
    A(i+(j-1)*m,i+(j-2)*m)= (2/k);
    A(i+(j-1)*m,i+(j-3)*m)= -(1/(2*k)); 
    b(i+(j-1)*m)=g2(x(i));
end

for j=1:n % left(power) and right boundary points
    i=1;
    A(i+(j-1)*m,i+(j-1)*m) = (-3/(2*h));
    A(i+(j-1)*m,(i+1)+(j-1)*m) = (2/h);
    A(i+(j-1)*m,(i+2)+(j-1)*m) = -(1/(2*h));
    if (j - 1)*k >= L_start && (j -1 )* k <= L_stop && L>0
        b(i+(j-1)*m)=P(y(j));
    else
        b(i+(j-1)*m)=g3(y(j));
    end
    i=m;
    A(i+(j-1)*m,i+(j-1)*m)=((-3/(2*h))+(H/K));
    A(i+(j-1)*m,(i-1)+(j-1)*m)=(2/h);
    A(i+(j-1)*m,(i-2)+(j-1)*m)=-(1/(2*h));
    if (j - 1)*k >= R_start && (j -1 )* k <= R_stop && L2>0
        b(i+(j-1)*m)=P2(y(j));
    else
        b(i+(j-1)*m)=g4(y(j));
    end
end

v=(A\b)+20; % solve for solution in v labeling
w=reshape(v(1:mn),m,n); %translate from v to w
if plot

mesh(x,y,w')
title("Heat distribution of a heatsink")
grid on
zlabel("T[°C]")
xlabel("x [cm]")
ylabel("y [cm]")
end