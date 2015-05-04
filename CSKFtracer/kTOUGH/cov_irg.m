function D=cov_irg(varargin)

% function to calculate covariance matrix of irregular grid
% Amalia Kokkinaki
% Sep 5, include anisotropy


% x is the VECTOR of all x's 
% y is the VECTOR of all y's
% e.g. for a 2x2 grid  x=[1 1 2 2] y=[1 2 1 2]

% if grid is regular this can be done by: 
% [x,y]=meshgrid(x,y);
% [xj,xl]=meshgrid(x(:),x(:));
% [yj,yl]=meshgrid(y(:),y(:));
% D=cov_irg(xl(:,1),yl(:,1));

if nargin == 4
    % 2D problem
    x=varargin{1};
    y=varargin{2};
    z=zeros(size(y));
    lx = varargin{3};
    ly = varargin{4};
    lz = 1;
elseif nargin == 6
    % 3D problem
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    lx = varargin{4};
    ly = varargin{5};
    lz = varargin{6};
end

m=length(x);
if (m-length(y)) || (m-length(z))
    disp('x, y and z need to be the same size')
end

D=zeros(m,m);
tic

for k=1:m
    for l=1:m
	rx = ((x(k)-x(l))/lx).^2;
	ry = ((y(k)-y(l))/ly).^2;
	rz = ((z(k)-z(l))/lz).^2;
	D(k,l) = sqrt(rx + ry + rz);
%       D(k,l)=sqrt((x(k)-x(l)).^2 + (y(k)-y(l)).^2 + (z(k)-z(l)).^2);
    end
end

toc
end

function D2=cov_reg(varargin)

% function to calculate covariance for regular grid using vectorized form
% to check how slower my code is

[x,y,z] = meshgrid(x,y,z);
[xj,xl]=meshgrid(x(:),x(:));
[yj,yl]=meshgrid(y(:),y(:));
[zj,zl]=meshgrid(z(:),z(:));
D2 = sqrt(((xj-xl)).^2+((yj-yl)).^2+((zj-zl)).^2);

end
