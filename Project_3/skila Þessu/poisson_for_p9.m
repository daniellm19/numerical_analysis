function w=poisson2(xl,xr,yb,yt,M,N,L_start,L_stop,p, K, H, plot, rect_size_points)
% Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps
% Output: matrix w holding solution values
% Example usage: w=poisson(0,1,1,2,4,4)

%Power Variables
delta = 0.1;
L = L_stop - L_start;

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

S = zeros(m,n);

% Create matrix for identifying forbidden nodes
if nargin < 13
    rect_size_points = 0
else
    rect_size_points = rect_size_points
    rect_points = zeros(rect_size_points,rect_size_points)
    for i=m:-1:(m-rect_size_points+1)
        for j=n:-1:(n-rect_size_points+1)
            S(i,j) = 1;
        end
    end
end

for i=2:m-1 % interior points
    for j=2:n-1
        loc = i+(j-1)*m;
        if S(i,j) == 1
            A(loc,loc)=0.000000001;
            continue     % Landed on forbidden node
        end
        if S(i+1,j) == 1 % Apply right if left to forbidden box
            S(i,j) = 3;
            A(loc,loc)=       (((2*h*H)/K)-3);
            A(loc,loc-1)=     4;
            A(loc,loc-2)=     -1;
            continue
        end
        if S(i,j+1) == 1 % Apply top if bottom to forbidden box
            S(i,j) = 4;
            A(loc,loc)=     (((2*k*H)/K)-3);
            A(loc,loc-m)=     4;
            A(loc,loc-2*m)=     -1; 
            continue
        end
        if S(i+1,j+1) == 1 % Apply top if bottom to forbidden box
            S(i,j) = 4;
            A(loc,loc)=     (((2*k*H)/K)-3);
            A(loc,loc-m)=     4;
            A(loc,loc-2*m)=     -1; 
            continue
        end

        S(i,j) = 6;
        A(loc,loc)=     -2*((1/h2)+(1/k2)+(H/(K*delta)));
        A(loc,loc+1)=   1/h2;
        A(loc,loc-1)=   1/h2;
        A(loc,loc-m)=   1/k2;
        A(loc,loc+m)=   1/k2;
    end
end
for i=2:m-1 % bottom and top boundary points
    j=1;    % Bottom
    loc = i+(j-1)*m;
    S(i,j) = 5;
    A(loc,loc)=     (((2*k*H)/K)-3);
    A(loc,loc+m) =        4;
    A(loc,loc+2*m) =    -1;

    j=n;    % Top
    loc = i+(j-1)*m;
    if S(i,j) == 1
        A(loc,loc)=0.000000001;
        continue     % Landed on forbidden node
    end
    if S(i+1,j) == 1 % Apply right if left to forbidden box
        S(i,j) = 3;
        A(loc,loc)=       (((2*h*H)/K)-3);
        A(loc,loc-1)=     4;
        A(loc,loc-2)=     -1;
        continue
    end
    S(i,j) = 4;
    A(loc,loc)=     (((2*k*H)/K)-3);
    A(loc,loc-m)=     4;
    A(loc,loc-2*m)=     -1; 
end

for j=1:n % left(power) and right boundary points
    i=1;    % left
    loc = i+(j-1)*m;
    S(i,j) = 2;
    A(loc,loc) =      -3;
    A(loc,loc+1) =    4;
    A(loc,loc+2) =    -1;
    b(loc) =          -((2*h*p)/(L*delta*K));
    i=m;    % right
    loc = i+(j-1)*m;
    if S(i,j) == 1
        A(loc,loc)=0.0000001;
        continue     % Landed on forbidden node
    end
%     if (j~=n)
%         if (S(i,j+1) == 1) % Apply top if bottom to forbidden box
%             S(i,j) = 4;
%             A(loc,loc)=     (((2*k*H)/K)-3);
%             A(loc,loc-m)=     4;
%             A(loc,loc-2*m)=     -1; 
%             continue
%         end
%     end
    S(i,j) = 3;
    A(loc,loc)=       (((2*h*H)/K)-3);
    A(loc,loc-1)=     4;
    A(loc,loc-2)=     -1;
end
A
S
v=(A\b)+20; % solve for solution in v labeling
w=reshape(v(1:mn),m,n); %translate from v to w
for i=m:-1:(m-rect_size_points+1)
    for j=n:-1:(n-rect_size_points+1)
        w(i,j) = NaN;
    end
end
if plot
mesh(x,y,w')
end