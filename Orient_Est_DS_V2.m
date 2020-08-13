function [azimuth,elevation,circularity] = Orient_Est_DS_V2(vol,W_est,W_circ,Thresh_circ,res)
% function determines the dominate 2D orientation of texture in the volume 
% im. The computation is performed in a local window of size W by W by W.

%%
% Define size and sigma for Gaussian filters
s = floor(W_est/2);      % Size of 3D Gaussian filters (filters are s by s by 2 voxels in size)
sigma = (s+2)/4;  % Sigma value of the 3D Gaussians (assuming symmetrical Gaussians)

% determine size of image:
[M,N,P] = size(vol);


%%
%--------------------------------------------------------------------------
% Create Gaussian Filters: 
%--------------------------------------------------------------------------
% Define 1D filters:
g=@(x)exp(-x.*x/sigma^2);     % 1D Gaussian function
G= g(-s:s);
Gd= (-s:s).*G;      % Derivative of a Guassian
Gd=Gd.'/norm((-s:s).*Gd,1);
G=G.'/norm(G,1);


%%
%--------------------------------------------------------------------------
% Calculate image gradients and determine per pixel estimate of orientation 
%--------------------------------------------------------------------------

% Compute gradients using filtering:
Gx = Quick3DFilter(vol,Gd,G,1);
Gy = Quick3DFilter(vol,Gd,G,2);
Gz = Quick3DFilter(vol,Gd,G,3);

% calculate magnitude and angles for gradient assuming it is a 3D vector:
G = single(sqrt(abs(Gx).^2 + abs(Gy).^2 + abs(Gz).^2));      % Magnitude
% Azimuth angle (in radians):
phi = single(atan2(Gx,Gy));                 % Note that x and y are swapped due to how Matlab defines coordinates
phi(phi<0) = phi(phi<0)+2*pi;
% Elevation angle (in radians):
theta = single(acos(Gz./G));

%%
%--------------------------------------------------------------------------
% Determine whether each of the planes is circular in local window:
%--------------------------------------------------------------------------
Gxy = gpuArray(Gx + 1j*Gy);
Gxz = gpuArray(Gx + 1j*Gz);
Gyz = gpuArray(Gy + 1j*Gz);

Rxy = local_average(abs(Gxy).^2, W_circ);
Rxz = local_average(abs(Gxz).^2, W_circ);
Ryz = local_average(abs(Gyz).^2, W_circ);

R1xy = local_average(Gxy.^2, W_circ);
R1xz = local_average(Gxz.^2, W_circ);
R1yz = local_average(Gyz.^2, W_circ);

Circularity_XY = abs(R1xy).^2./Rxy.^2;
Circularity_XZ = abs(R1xz).^2./Rxz.^2;
Circularity_YZ = abs(R1yz).^2./Ryz.^2;

circularity = zeros(M,N,P);
circularity(Circularity_XY>Thresh_circ | Circularity_XZ>Thresh_circ | Circularity_YZ>Thresh_circ) = 1;

%%
%--------------------------------------------------------------------------
% Determine dominant orientations in local window: 
%--------------------------------------------------------------------------

% determine a vector of possible orientations:
angle_vector2 = single([0,res/2,res:res:(180-res),(180-res/2)]);
N_angle2 = numel(angle_vector2);
angle_vector1 = single([0,res/2,res:res:(180-res),(180-res/2)]);
N_angle1 = numel(angle_vector1);


%% 2D Search
min_Values = ones(M,N,P, 'single', 'gpuArray')*100000;
dominant_theta = -ones(M,N,P, 'single', 'gpuArray');
dominant_phi = -ones(M,N,P, 'single', 'gpuArray');
% for loop to calculate the individual images for each value of theta:
for index1 = 1:N_angle1
    fprintf('Theta Angle: %d\n',index1);
    theta_est = angle_vector1(index1)*pi/180;
    for index2 = 1:N_angle2
        % obtain local value of phi and convert to radians:
        phi_est = angle_vector2(index2)*pi/180;
    
        % calculate test:
        test_Values = G.*abs(cos(theta_est).*cos(theta) + sin(theta_est).*sin(theta).*cos((phi_est-phi)));
        testGPU = gpuArray(test_Values);

        % perform local average:
        testGPU = local_average(testGPU, W_est);

        % find the values to replace
        [min_Values,dominant_theta,dominant_phi] = arrayfun(@replace_Values,...
            min_Values,testGPU,dominant_theta,dominant_phi,theta_est,phi_est);
    end

end

% Gather results from GPU:
azimuth = gather(dominant_phi);
elevation = gather(dominant_theta);
end

%% Embedded functions 
function im=local_average(im,W)
%% Performs averaging over local region
im = movmean(im, W, 1);
im = movmean(shiftdim(im,1), W, 1);
im = movmean(shiftdim(im,1), W, 1);
im = shiftdim(im,1);
end

function [orig,theta,phi] = replace_Values(orig, test, theta, phi, theta_est, phi_est)
%% Find values where the test data is less than the minimum data
if test < orig
    orig = test;
    phi = phi_est;
    theta = theta_est;
end
end

function G_out = Quick3DFilter(im,Gd,G,index)
%% Efficient 3D convolution for matlab
switch index
    case 1
        % hx filtering:
        G_out = imfilter(im, G,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), Gd,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), G,'symmetric'); 
        G_out = shiftdim(G_out,1);
    case 2
        % hy filtering
        G_out = imfilter(im, Gd,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), G,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), G,'symmetric'); 
        G_out = shiftdim(G_out,1);
    case 3
        % hz filtering
        G_out = imfilter(im, G,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), G,'symmetric'); 
        G_out = imfilter(shiftdim(G_out,1), Gd,'symmetric'); 
        G_out = shiftdim(G_out,1);
end
end