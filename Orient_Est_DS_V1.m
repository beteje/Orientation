function Angle = Orient_Est_DS_V1(im,W,num_scales,trig_average)

if nargin == 2
    num_scales = 1;
    trig_average = 0;
elseif nargin == 3
    trig_average = 0;
end

num_scales = 1;
trig_average = 0;

%%
% Define size and sigma for Gaussian filters
s = 3;      % Size of 3D Gaussian filters (filters are s by s by 2 voxels in size)
sigma = (s+2)/4;  % Sigma value of the 3D Gaussians (assuming symmetrical Gaussians)

% determine size of image:
[M,N,P] = size(im);


%%
%--------------------------------------------------------------------------
% Create Gaussian Filters: 
%--------------------------------------------------------------------------
% disp('Calculating Orientation of Stress Fibres....');
% Create i and j coordinates for the filters:
% [i,j,k] = meshgrid(-s:s,-s:s,-s:s);
% % Create filters:
% hx = (2.*i./sigma^2).*exp( -(i.^2 + j.^2 + k.^2)./sigma^2);    % filter in x
% hy = (2.*j./sigma^2).*exp( -(i.^2 + j.^2 + k.^2)./sigma^2);    % filter in y
% hz = (2.*k./sigma^2).*exp( -(i.^2 + j.^2 + k.^2)./sigma^2);    % filter in z
% 
% hx=hx/sum(sum(sum(i.*hx)));
% hy=hy/sum(sum(sum(j.*hy)));
% hz=hz/sum(sum(sum(k.*hz)));

% Define 1D filters:
g=@(x)exp(-x.*x/sigma^2);     % 1D Gaussian function
G= g(-s:s);
Gd= (-s:s).*G;      % Derivative of a Guassian
Gd=Gd/norm((-s:s).*Gd,1);
G=G/norm(G,1);


%%
%--------------------------------------------------------------------------
% Calculate image gradients and determine per pixel estimate of orientation 
%--------------------------------------------------------------------------

% Compute gradients using filtering:

% Gx = imfilter(im,hx,'symmetric');   % Gradient in x
% Gy = imfilter(im,hy,'symmetric');   % Gradient in y
% Gz = imfilter(im,hz,'symmetric');   % Gradient in z

% More efficient filtering:
Gx = Quick3DFilter(im,Gd,G,1);
Gy = Quick3DFilter(im,Gd,G,2);
Gz = Quick3DFilter(im,Gd,G,3);

% Determine normal vector field assuming Girdle distribution:
[Px,Py,Pz] = DeterminePolarAxis(Gx,Gy,Gz,W);

% calculate magnitude and angles for gradient assuming it is a 3D vector:
G = sqrt(abs(Gx).^2 + abs(Gy).^2 + abs(Gz).^2);      % Magnitude
% Azimuth angle (in radians):
phi = atan2(Px,Py);                 % Note that x and y are swapped due to how Matlab defines coordinates
% Elevation angle (in radians):
% theta = acos(Pz);
theta = abs(atan(sqrt(Px.^2 + Py.^2)./Pz));


%%
%--------------------------------------------------------------------------
% Determine dominant orientations in local window: 
%--------------------------------------------------------------------------

% determine a vector of possible orientations:
res = 2;
angle_vector1 = 0:res:(360-res);
N_angle1 = numel(angle_vector1);
angle_vector2 = 0:res:(90-res);
% angle_vector2 = -90:res:(90-res);
N_angle2 = numel(angle_vector2);
% 
% angle_vector1 = 0:res:(360-res);
% N_angle1 = numel(angle_vector1);
% angle_vector2 = 0:res:(90-res);
% N_angle2 = numel(angle_vector2);

kappa = 2;

% %% 2D Search
% % Initialise 5D Matrix to hold data
% A = zeros(M,N,P,N_angle,N_angle);
% % for loop to calculate the individual images for each value of theta:
% for index1 = 1:N_angle
%     theta_est = angle_vector(index1)*pi/180;
%     for index2 = 1:N_angle
%            
%         % obtain local value of theta and conver to radians:
%         phi_est = angle_vector(index2)*pi/180;
%     
%         % calculate test:
%         holder = 2*kappa.*G.*(1 - cos(theta_est).*cos(theta) - sin(theta_est).*sin(theta).*cos(2.*(phi_est-phi)));
%     
%         % perform local summation:
%         A(:,:,:,index1,index2) = local_sum(holder,W);
%     end
% 
% end
% 
% [~,i2] = min(reshape(A,M,N,P,N_angle*N_angle),[],4);
% [i1,i2] = ind2sub([N_angle,N_angle],i2(:));

%% 2 1D Searches
% Initialise 4D Matrix to hold data
A = zeros(M,N,P,N_angle1);

% Determine phi based on theta_est=90 degrees:
for index1 = 1:N_angle1
    % obtain local value of theta and conver to radians:
    phi_est = angle_vector1(index1)*pi/180;
    % calculate test:
%     holder = 2*kappa.*G.*(1 - sin(theta).*cos(2.*(phi_est-phi)));
    holder = 2*kappa.*G.*(1 - sin(theta).*cos((phi_est-phi)));
    % perform local summation:
    A(:,:,:,index1,1) = local_sum(holder,W);
end
[~,i2] = min(A,[],4);
Dom_angle = angle_vector1(i2);

% Interp peak:
Delta=subSample(A,i2,res);
Dom_angle = Dom_angle+reshape(Delta,M,N,P);

A = zeros(M,N,P,N_angle2);
% phi_est = Dom_angle.*pi./180;
phi_est = pi/5;
for index1 = 1:N_angle2
    % obtain local value of theta and conver to radians:
    theta_est = angle_vector2(index1)*pi/180;
    % calculate test:
    holder = 2*kappa.*G.*(1 - cos(theta_est).*cos(theta) - sin(theta_est).*sin(theta).*cos((phi_est-phi)));
%     holder = 2*kappa.*G.*(1 - cos(theta_est).*cos(theta) - sin(theta_est).*sin(theta));
    % perform local summation:
    A(:,:,:,index1) = local_sum(holder,W);
end
[~,i2] = min(A,[],4);
Dom_angle(:,:,:,2) = angle_vector2(i2);

% Interp peak:
Delta=subSample(A,i2,res);
Dom_angle(:,:,:,2) = Dom_angle(:,:,:,2)+reshape(Delta,M,N,P);

% Convert to radian:
Angle = Dom_angle.*pi/180;

if trig_average == 1
    Angle = angle(local_sum(exp(1i.*2.*Angle),W))./2;
end


end

%% Embedded functions 
function J=local_sum(im,W)
%% Performs summation over local region
% J=imfilter(im,ones(W,W,W),'symmetric');
F=ones([1 W]);
FA = reshape(F,[length(F) 1 1]);
J = imfilter(im, FA,'symmetric');
J = imfilter(shiftdim(J,1), FA,'symmetric');
J = imfilter(shiftdim(J,1), FA,'symmetric');
J = shiftdim(J,1)/(W*W*W);
end

function Delta=subSample(A,index2,res)
%% Determine subsample accuracy of a peak
% Size
[M,N,P,Q] = size(A);

% indexing coordinates:
[Y,X,Z] = meshgrid(1:M,1:N,1:P);

% Interp peak:
index1 = index2-1;
index1(index1==0) = Q;
R1 = A(sub2ind([M,N,P,Q],X(:),Y(:),Z(:),index1(:)));
R2 = A(sub2ind([M,N,P,Q],X(:),Y(:),Z(:),index2(:)));
index3 = index2+1;
index3(index3==Q+1) = 1;
R3 = A(sub2ind([M,N,P,Q],X(:),Y(:),Z(:),index3(:)));

% Gaussian Interp:
Delta=((log(R3(:))-log(R1(:)))./(4*log(R2(:))-2*log(R1(:))-2*log(R3(:)))).*res;
% % parabola Interp:
% Delta=((R3(:)-R1(:))./(2*(2*R2(:)-R1(:)-R3(:)))).*res;

Delta(isnan(Delta)) = 0;
end

function G_out = Quick3DFilter(im,Gd,G,index)
%% Efficient 3D convolution for matlab
switch index
    case 1
        % hx filtering:
        G_out = imfilter(im, reshape(G,[length(Gd) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(Gd,[length(G) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
        G_out = shiftdim(G_out,1);
    case 2
        % hy filtering
        G_out = imfilter(im, reshape(Gd,[length(Gd) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
        G_out = shiftdim(G_out,1);
    case 3
        % hz filtering
        G_out = imfilter(im, reshape(G,[length(Gd) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), reshape(Gd,[length(G) 1 1]),'symmetric'); 
        G_out = shiftdim(G_out,1);
end
end

function e1 = SmallestEig(a11, a12, a13, a22, a23, a33)
%% function determines the smallest eigenvalues for the symmetric 3 by 3 
% matrix defined by a11, a12, a13, a22, a23, a33

% determine coefficients for x^3 + bx^2 + cx + d= 0
b = -(a11+a22+a33);
q = -b./3;      % trace of the matrix
h1 = a12.^2;
h2 = a13.^2;
h3 = a23.^2;
p1 = h1+h2+h3;
c = -(p1 - a22.*a33 - a11.*a22 - a11.*a33);
d = -(a11.*a22.*a33 - a11.*h3 - a22.*h2 - a33.*h1 + 2.*a12.*a13.*a23);
clear('h1','h2','h3');


p2 = (b.^2).*2./3 - 2.*c;
p = sqrt(p2./6);
r = -(q.^2.*(q + b) + c.*q + d)./(p.^3.*2);

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
phi = acos(r)./3;
phi(r<=-1) = pi./3;
phi(r>=1) = 0;

% smallest eigenvalue:
e1 = q + 2.*p.*cos(phi + (2*pi/3));
end

function [Px_out,Py_out,Pz_out] = DeterminePolarAxis(Gx,Gy,Gz,W)
%% Determine polar axis for gradient vector field by assuming the local 
% field follows a Girdle distribution. Thus the eigenvector corresponding
% to the smallest eigenvalue of the scatter matrix defines the polar axis
% (should be perpendicular to the local gradient vectors).

% determine size of Gx
[s1,s2,s3] = size(Gx);

% Determine scatter matrices (i.e. structure tensor):
a11 = local_sum(Gx.^2,W);
a22 = local_sum(Gy.^2,W);
a33 = local_sum(Gz.^2,W);
a12 = local_sum(Gx.*Gy,W);
a13 = local_sum(Gx.*Gz,W);
a23 = local_sum(Gy.*Gz,W);

% Determine smallest eigenvalue:
e1 = SmallestEig(a11(:), a12(:), a13(:), a22(:), a23(:), a33(:));

% Determine eigenvector assuming x component is 1:
A=zeros(s1*s2*s3,2,2);            
b=zeros(s1*s2*s3,2);          

% % Set up 2 by 2 matrix per voxel:
% A(:,1,1) = a12(:);
% A(:,1,2) = a13(:);
% A(:,2,1) = a22(:)-e1(:);
% A(:,2,2) = a23(:);
% % Set up 2 by 1 vector per voxel:
% b(:,1) = -a11(:)+e1(:);
% b(:,2) = -a12(:);
% 
% % Gauss elimination to determine eigenvector:
% Nb = 3;
% coeffs=zeros(s1*s2*s3,Nb-1);            %single/double
% for k=1:(Nb-1)
%     for l=(k+1):(Nb-1)
%         c=A(:,l,k)./A(:,k,k);
%         for m=(k+1):(Nb-1)
%             A(:,l,m)=A(:,l,m)-c.*A(:,k,m);
%         end
%         A(:,l,k)=0;
%         b(:,l)=b(:,l)-c.*b(:,k);
%     end
% end
% for k=(Nb-1):-1:1
%     coeffs(:,k)=b(:,k);
%     for m=(k+1):(Nb-1)
%         coeffs(:,k)=coeffs(:,k)-A(:,k,m).*coeffs(:,m);
%     end
%     coeffs(:,k)=coeffs(:,k)./A(:,k,k);
% end
% 
% % Determine length of vectors
% r = sqrt(1 + coeffs(:,1).^2 + coeffs(:,2).^2);
% 
% % Convert vectors to unit length and reshape:
% Px_out = reshape(1./r,s1,s2,s3);
% Py_out = reshape(coeffs(:,1)./r,s1,s2,s3);
% Pz_out = reshape(coeffs(:,2)./r,s1,s2,s3);


% Set up 2 by 2 matrix per voxel:
A(:,1,1) = a11(:)-e1(:);
A(:,1,2) = a12(:);
A(:,2,1) = a12(:);
A(:,2,2) = a22(:)-e1(:);
% Set up 2 by 1 vector per voxel:
b(:,1) = -a13(:);
b(:,2) = -a23(:);

% Gauss elimination to determine eigenvector:
Nb = 3;
coeffs=zeros(s1*s2*s3,Nb-1);            %single/double
for k=1:(Nb-1)
    for l=(k+1):(Nb-1)
        c=A(:,l,k)./A(:,k,k);
        for m=(k+1):(Nb-1)
            A(:,l,m)=A(:,l,m)-c.*A(:,k,m);
        end
        A(:,l,k)=0;
        b(:,l)=b(:,l)-c.*b(:,k);
    end
end
for k=(Nb-1):-1:1
    coeffs(:,k)=b(:,k);
    for m=(k+1):(Nb-1)
        coeffs(:,k)=coeffs(:,k)-A(:,k,m).*coeffs(:,m);
    end
    coeffs(:,k)=coeffs(:,k)./A(:,k,k);
end

% Determine length of vectors
r = sqrt(1 + coeffs(:,1).^2 + coeffs(:,2).^2);

% Convert vectors to unit length and reshape:
Pz_out = reshape(1./r,s1,s2,s3);
Px_out = reshape(coeffs(:,1)./r,s1,s2,s3);
Py_out = reshape(coeffs(:,2)./r,s1,s2,s3);


% set null cases to zero:
index = isnan(e1);
Px_out(index) = 0;
Py_out(index) = 0;
Pz_out(index) = 0;

end

function [Px,Py,Pz] = PolarAxis_test(Gx,Gy,Gz,W)
%% Similar to DeterminePolarAxis but computes eigenvectors and eigenvalues
% directly per voxel using a for loop

[s1,s2,s3] = size(Gx);
% Determine scatter matrices (i.e. structure tensor):
a11 = local_sum(Gx.^2,W);
a22 = local_sum(Gy.^2,W);
a33 = local_sum(Gz.^2,W);
a12 = local_sum(Gx.*Gy,W);
a13 = local_sum(Gx.*Gz,W);
a23 = local_sum(Gy.*Gz,W);

A(:,1,1) = a11(:);
A(:,1,2) = a12(:);
A(:,1,3) = a13(:);
A(:,2,1) = a12(:);
A(:,2,2) = a22(:);
A(:,2,3) = a23(:);
A(:,3,1) = a13(:);
A(:,3,2) = a23(:);
A(:,3,3) = a33(:);

coeffs = zeros(s1*s2*s3,3);
for index = 1:s1*s2*s3
    [V,D] = eig(squeeze(A(index,:,:)));
    [~,i2] = min(diag(D));
    
    coeffs(index,:) = V(:,i2).';
end

Px = reshape(coeffs(:,1),s1,s2,s3);
Py = reshape(coeffs(:,2),s1,s2,s3);
Pz = reshape(coeffs(:,3),s1,s2,s3);
end