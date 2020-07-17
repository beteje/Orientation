function Angle = Orient_Est_DS_V2(im,W)
% function determines the dominate 2D orientation of texture in the volume 
% im. The computation is performed in a local window of size W by W by W.

%%
% Define size and sigma for Gaussian filters
s = 6;      % Size of 3D Gaussian filters (filters are s by s by 2 voxels in size)
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


% calculate magnitude and angles for gradient assuming it is a 3D vector:
G = sqrt(abs(Gx).^2 + abs(Gy).^2 + abs(Gz).^2);      % Magnitude
% Azimuth angle (in radians):
phi = atan2(Gx,Gy);                 % Note that x and y are swapped due to how Matlab defines coordinates
phi(phi<0) = phi(phi<0)+2*pi;
% Elevation angle (in radians):
theta = acos(Gz./G);


%%
%--------------------------------------------------------------------------
% Determine dominant orientations in local window: 
%--------------------------------------------------------------------------

% determine a vector of possible orientations:
res = 1;
angle_vector2 = [0,res/2,res:res:(180-res),(180-res/2)];
N_angle2 = numel(angle_vector2);
angle_vector1 = [0,res/2,res:res:(180-res),(180-res/2)];
N_angle1 = numel(angle_vector1);


%% 2D Search
% % Initialise 5D Matrix to hold data
A = zeros(M,N,P,N_angle1,N_angle2);
% for loop to calculate the individual images for each value of theta:
for index1 = 1:N_angle1
    theta_est = angle_vector1(index1)*pi/180;
    for index2 = 1:N_angle2
           
        % obtain local value of theta and conver to radians:
        phi_est = angle_vector2(index2)*pi/180;
    
        % calculate test:
        holder = G.*abs(cos(theta_est).*cos(theta) + sin(theta_est).*sin(theta).*cos((phi_est-phi)));
    
        % perform local summation:
        A(:,:,:,index1,index2) = local_sum(holder,W);
    end

end

% find the minimum solution:
[~,i2] = min(reshape(A,M,N,P,N_angle1*N_angle2),[],4);
[i1,i2] = ind2sub([N_angle1,N_angle2],i2(:));
Dom_angle_phi = reshape(angle_vector2(i2),M,N,P);           % phi
Dom_angle_theta = reshape(angle_vector1(i1),M,N,P);         % theta

% Determine sub pixel accuracy:
[Delta_theta,Delta_phi]=subSample2D(A,i1,i2,res,res);
Dom_angle_phi = Dom_angle_phi + reshape(Delta_phi,M,N,P);
Dom_angle_theta = Dom_angle_theta + reshape(Delta_theta,M,N,P);

% Deal with any cases where phi<0:
index = Dom_angle_phi<0;
Dom_angle_phi(index) = Dom_angle_phi(index)+180;
Dom_angle_theta(index) = 180-Dom_angle_theta(index);

% Convert to radian:
Angle = Dom_angle_phi.*pi/180;
Angle(:,:,:,2) = Dom_angle_theta.*pi/180;


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

function [Delta_x,Delta_y]=subSample2D(Data,index_x,index_y,res_x,res_y)
%% Determine subsample accuracy of a peak based on 2D fitting using the
% form y = Ax^2 + By^2 + Cx + Dy + Exy + F

% Size of A
[M,N,P,Q,R] = size(Data);

% indexing coordinates:
[Y,X,Z] = meshgrid(1:M,1:N,1:P);

% Define x indices:
index_x11 = index_x-1;
index_x11(index_x11==0) = Q;
index_x12 = index_x11;
index_x13 = index_x11;
index_x21 = index_x;
index_x23 = index_x;
index_x31 = index_x+1;
index_x31(index_x31==Q+1) = 1;
index_x32 = index_x31;
index_x33 = index_x31;

% Define y indices:
index_y1 = index_y-1;
holder1 = logical(index_y1==0);
index_y1(holder1==1) = R;
index_x11(holder1==1) = Q-index_x11(holder1==1)+2;
index_x21(holder1==1) = Q-index_x21(holder1==1)+2;
index_x31(holder1==1) = Q-index_x31(holder1==1)+2;
index_x11(index_x11==Q+1) = 1;
index_x21(index_x21==Q+1) = 1;
index_x31(index_x31==Q+1) = 1;

index_y3 = index_y+1;
holder1 = logical(index_y3==(R+1));
index_y3(holder1==1) = 1;
index_x13(holder1==1) = Q-index_x13(holder1==1)+2;
index_x23(holder1==1) = Q-index_x23(holder1==1)+2;
index_x33(holder1==1) = Q-index_x33(holder1==1)+2;
index_x13(index_x13==Q+1) = 1;
index_x23(index_x23==Q+1) = 1;
index_x33(index_x33==Q+1) = 1;

% Obtain coefficients:
F = Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x(:),index_y(:)));

% determine A:
A = 0.5.*(Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x32(:),index_y(:))) + ...
    Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x12(:),index_y(:)))) - F;

% determine B:
B = 0.5.*(Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x23(:),index_y3(:))) + ...
    Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x21(:),index_y1(:)))) - F;

% determine C:
Z1 = Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x33(:),index_y3(:)));
C = 0.5.*(Z1 + Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x31(:),index_y1(:))))...
    -A-B-F;

% determine D:
D = 0.5.*(Z1 + Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x13(:),index_y3(:))))...
    -A-B-F;

% Determine E:
E = Data(sub2ind([M,N,P,Q,R],X(:),Y(:),Z(:),index_x11(:),index_y1(:))) - A- B + C + D - F;

% determine shifts:
Denom = 4.*A.*B - E.^2;

% x:
Delta_x = -(2.*B.*C-D.*E)./Denom.*res_x;
Delta_x(isnan(Delta_x)) = 0;
Delta_x(abs(Delta_x)>res_x)=0;

% y:
Delta_y = -(2.*A.*D-C.*E)./Denom.*res_y;
Delta_y(isnan(Delta_y)) = 0;
Delta_y(abs(Delta_y)>res_y)=0;
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