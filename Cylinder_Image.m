
M = 151;
N = 151;
P = 151;
M1 = floor(M/2);
N1 = floor(N/2);
P1 = floor(P/2);

% azimuth and elevation of cylinder:
azimuth1 = 20;
elevation1 = 40;

% function to indecies of voxels in the line for a given x0, y0, azimuth
% and elevation
Line_fn = @(x0,y0, elevation1,azimuth1) sub2ind([M,N,P],round(M1+((-(M1-1):(M1-1)) - x0).*sin(elevation1*pi./180)*cos(azimuth1.*pi./180)),round(N1+((-(N1-1):(N1-1)) - y0).*sin(elevation1*pi./180)*sin(azimuth1.*pi./180)), round(P1+(-(P1-1):(P1-1)).*cos(elevation1*pi./180)));

% construct image
A = zeros(151,151,151);
A(Line_fn(0,0,elevation1,azimuth1)) = 1;
A(Line_fn(1,0,elevation1,azimuth1)) = 1;
A(Line_fn(0,1,elevation1,azimuth1)) = 1;
A(Line_fn(-1,0,elevation1,azimuth1)) = 1;
A(Line_fn(0,-1,elevation1,azimuth1)) = 1;

% % include these to turn cylinder into a cubiod
% A(Line_fn(1,1,elevation1,azimuth1)) = 1;
% A(Line_fn(-1,1,elevation1,azimuth1)) = 1;
% A(Line_fn(1,-1,elevation1,azimuth1)) = 1;
% A(Line_fn(-1,-1,elevation1,azimuth1)) = 1;

% Blur image with Gaussian filter:
s = 5; 
sigma = 2;
[i,j,k] = meshgrid(-s:s,-s:s,-s:s);
A_im = imfilter(double(A),exp( -(i.^2 + j.^2 + k.^2)./sigma^2));