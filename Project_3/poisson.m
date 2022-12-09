function w=poisson(xl,xr,yb,yt,M,N,L,plot)
% Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps
% Output: matrix w holding solution values
% Example usage: w=poisson(0,1,1,2,4,4)

%Power Variables
p = 5;
K = 1.68;
H = 0.005;
delta = 0.1;

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
        A(i+(j-1)*m,i+(j-1)*m)=     -2*((1/h2)+(1/k2)+(H/(K*delta)));
        A(i+(j-1)*m,i+1+(j-1)*m)=   1/h2;
        A(i+(j-1)*m,(i-1)+(j-1)*m)= 1/h2;
        A(i+(j-1)*m,i+(j-2)*m)=     1/k2;
        A(i+(j-1)*m,i+j*m)=         1/k2;
    end
end
for i=2:m-1 % bottom and top boundary points
    j=1;
    A(i+(j-1)*m,i+(j-1)*m)=     (((2*k*H)/K)-3);
    A(i+(j-1)*m,i+j*m) =        4;
    A(i+(j-1)*m,i+(j+1)*m) =    -1;

    j=n;
    A(i+(j-1)*m,i+(j-1)*m)=     (((2*k*H)/K)-3);
    A(i+(j-1)*m,i+(j-2)*m)=     4;
    A(i+(j-1)*m,i+(j-3)*m)=     -1; 
end

for j=1:n % left(power) and right boundary points
    i=1;
    A(i+(j-1)*m,i+(j-1)*m) =        -3;
    A(i+(j-1)*m,(i+1)+(j-1)*m) =    4;
    A(i+(j-1)*m,(i+2)+(j-1)*m) =    -1;
    b(i+(j-1)*m) =                  -((2*h*p)/(L*delta*K));

    i=m;
    A(i+(j-1)*m,i+(j-1)*m)=         (((2*h*H)/K)-3);
    A(i+(j-1)*m,(i-1)+(j-1)*m)=     4;
    A(i+(j-1)*m,(i-2)+(j-1)*m)=     -1;
end

v=(A\b)+20; % solve for solution in v labeling
w=reshape(v(1:mn),m,n); %translate from v to w
if plot
mesh(x,y,w')
end