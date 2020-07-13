M = 51;
N = 51;
P = 51;
M1 = floor(M/2);
N1 = floor(N/2);
P1 = floor(P/2);

% function to indecies of voxels in the line for a given x0, y0, azimuth
% and elevation
Line_fn = @(x0,y0, elevation1,azimuth1) sub2ind([M,N,P],round(M1+((-(M1-2):(M1-2)) - x0).*sin(elevation1*pi./180)*cos(azimuth1.*pi./180)),round(N1+((-(N1-2):(N1-2)) - y0).*sin(elevation1*pi./180)*sin(azimuth1.*pi./180)), round(P1+(-(P1-2):(P1-2)).*cos(elevation1*pi./180)));

% Blur image with Gaussian filter:
s = 5; 
sigma = 2;
[i,j,k] = meshgrid(-s:s,-s:s,-s:s);

PSNR = 20;

azimuth_array = rand(1).*180.*ones(1,500);
elevation_array = (rand(1,10)./1 + 0.00).*180;

% azimuth_array = 178 + rand;
% elevation_array = 50+rand;

s_array = [13];

phi_estZ = zeros(length(s_array),length(azimuth_array),length(elevation_array));
theta_estZ = zeros(length(s_array),length(azimuth_array),length(elevation_array));

phi_estY = zeros(length(s_array),length(azimuth_array),length(elevation_array));
theta_estY = zeros(length(s_array),length(azimuth_array),length(elevation_array));

phi_estT = zeros(length(s_array),length(azimuth_array),length(elevation_array));
theta_estT = zeros(length(s_array),length(azimuth_array),length(elevation_array));


for index_outer = 1:length(s_array)
for index_a = 1:length(azimuth_array)
    azimuth1 = azimuth_array(index_a);
    for index_e = 1:length(elevation_array)
        disp(['index_a = ', num2str(index_a/length(azimuth_array)*100,'%.0f')...
            ,'%, index_e = ', num2str(index_e/length(elevation_array)*100,'%.0f'),'%']); 
        elevation1 = elevation_array(index_e);
        
        A = imagebuild(azimuth1,elevation1,Line_fn,M,N,P);
%         A_im = imfilter(double(A),exp( -(i.^2 + j.^2 + k.^2)./sigma^2));
        A_im = QuickFilter(double(A));
        en = randn(size(A_im));
        sigma_noise = sqrt(max(abs(A_im(:)).^2)*10^(-PSNR/10)/mean(abs(en(:)).^2));
        A_im = A_im + en.*sigma_noise;
        index = logical(A>0);
        
%         Angle = Orient_Est_DS(A_im,15).*180./pi;
        
        %%%%%%%%%%%%%%
        %testing:
%         s = s_array(index_outer);
        s = 6;
        sigma = (s+2)/4;
        g=@(x)exp(-x.*x/sigma^2);     % 1D Gaussian function
        G= g(-s:s);
        Gd= (-s:s).*G;      % Derivative of a Guassian
        Gd=Gd/norm((-s:s).*Gd,1);
        G=G/norm(G,1);
        Gx = Quick3DFilter2(A_im,Gd,G,1);
        Gy = Quick3DFilter2(A_im,Gd,G,2);
        Gz = Quick3DFilter2(A_im,Gd,G,3);
        
        % calculate magnitude and angles for gradient assuming it is a 3D vector:
        G = sqrt(abs(Gx).^2 + abs(Gy).^2 + abs(Gz).^2);      % Magnitude
        % Azimuth angle (in radians):
        phi = atan2(Gx,Gy);                 % Note that x and y are swapped due to how Matlab defines coordinates
        phi(phi<0) = phi(phi<0)+2*pi;
        % Elevation angle (in radians):
        theta = acos(Gz./G);
        
        res1 = 0.5;
        res2 = 0.5;
        angle_vector2 = [0,res2/2,res2:res2:(180-res2),(180-res2/2)];
        N_angle2 = numel(angle_vector2);
        angle_vector1 = [0,res1/2,res1:res1:(180-res1),(180-res1/2)];
        % angle_vector2 = -90:res:(90-res);
        N_angle1 = numel(angle_vector1);
        
        W = s_array(index_outer);
        
        phi1 = phi((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        theta1 = theta((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        G1 = G((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        
        Gx1 = Gx((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        Gy1 = Gy((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        Gz1 = Gz((26-W):(26+W),(26-W):(26+W),(26-W):(26+W));
        
        
        % 2D Search
        % Initialise 5D Matrix to hold data
        Z = zeros(N_angle1,N_angle2);
%         Y = zeros(N_angle1,N_angle2);
        % for loop to calculate the individual images for each value of theta:
        for index1 = 1:N_angle1
            theta_est = angle_vector1(index1)*pi/180;
            for index2 = 1:N_angle2
                   
                % obtain local value of theta and conver to radians:
                phi_est = angle_vector2(index2)*pi/180;
            
                % calculate test:
                holder = G1.*abs(cos(theta1).*cos(theta_est) + sin(theta1).*sin(theta_est).*cos((phi1-phi_est)));
            
                % perform local summation:
                Z(index1,index2) = sum(holder(:))./((2*W+1)^3);
                
%                 holder = G1.*sqrt(1 - (cos(theta1).*cos(theta_est) + sin(theta1).*sin(theta_est).*cos((phi1-phi_est))).^2);
            
                % perform local summation:
%                 Y(index1,index2) = sum(holder(:))./W^3;
            end
        
        end
        
        [~,i2] = min(Z(:));
        [i1,i2] = ind2sub([N_angle1,N_angle2],i2);
        
        Z1(1,1,1,:,:) = Z;
        [Delta_theta,Delta_phi]=subSample2D(Z1,i1,i2,res1,res2);
        Angle_A1 = [angle_vector2(i2),angle_vector1(i1)]+[Delta_phi,Delta_theta];
        
        if Angle_A1(1)<0
            Angle_A1(1) = Angle_A1(1)+180;
            Angle_A1(2) = 180-Angle_A1(2);
        end
        
        [~,i2] = max(Y(:));
        [i1,i2] = ind2sub([N_angle1,N_angle2],i2);
        
%         Y1(1,1,1,:,:) = Y;
%         [Delta_theta1,Delta_phi1]=subSample2D(Y1,i1,i2,res1,res2);
%         Angle_B1 = [angle_vector2(i2),angle_vector1(i1)]+[Delta_phi1,Delta_theta1];
%         
%         if Angle_B1(1)<0
%             Angle_B1(1) = Angle_B1(1)+180;
%             Angle_B1(2) = 180-Angle_B1(2);
%         end
        
        T = [Gx1(:).';Gy1(:).';Gz1(:).']*[Gx1(:).';Gy1(:).';Gz1(:).']'./W^6;
        [V,D] = eigs(T,1,'smallestabs');
        Angle_T = [atan2(V(1),V(2))*180./pi, acos(V(3)).*180./pi];
        
        if Angle_T(1)<0
            Angle_T(1) = Angle_T(1)+180;
            Angle_T(2) = 180-Angle_T(2);
        end
                
        
        phi_estZ(index_outer,index_a,index_e) = Angle_A1(1);
        theta_estZ(index_outer,index_a,index_e) = Angle_A1(2);
%         phi_estY(index_outer,index_a,index_e) = Angle_B1(1);
%         theta_estY(index_outer,index_a,index_e) = Angle_B1(2);
        phi_estT(index_outer,index_a,index_e) = Angle_T(1);
        theta_estT(index_outer,index_a,index_e) = Angle_T(2);
        
    end
end
end

return

function A = imagebuild(az,el,line_fn,M,N,P)
    A = zeros(M,N,P); 
    index = line_fn(0,0,el,az);
    index(index>M*N*P) = M*N*P;
    index(index<1) = 1;
    A(index) = 1;
    
    index = line_fn(1,0,el,az);
    index(index>M*N*P) = M*N*P;
    index(index<1) = 1;
    A(index) = 1;
    
    index = line_fn(0,1,el,az);
    index(index>M*N*P) = M*N*P;
    index(index<1) = 1;
    A(index) = 1;
    
    index = line_fn(-1,0,el,az);
    index(index>M*N*P) = M*N*P;
    index(index<1) = 1;
    A(index) = 1;
    
    index = line_fn(0,-1,el,az);
    index(index>M*N*P) = M*N*P;
    index(index<1) = 1;
    A(index) = 1;
    
%     A(line_fn(1,1,el,az)) = 1;
%     A(line_fn(-1,1,el,az)) = 1;
%     A(line_fn(1,-1,el,az)) = 1;
%     A(line_fn(-1,-1,el,az)) = 1;
%     A(line_fn(-2,-2,el,az)) = 1;
%     A(line_fn(-2,2,el,az)) = 1;
%     A(line_fn(2,2,el,az)) = 1;
%     A(line_fn(2,-2,el,az)) = 1;
end

function G_out=QuickFilter(im)
%% Quick filter
s = 5;
sigma = 2;
g=@(x)exp(-x.*x/sigma^2);     % 1D Gaussian function
G= g(-s:s);
G=G/norm(G,1);
G_out = imfilter(im, reshape(G,[length(G) 1 1]),'symmetric'); 
G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
G_out = imfilter(shiftdim(G_out,1), reshape(G,[length(G) 1 1]),'symmetric'); 
G_out = shiftdim(G_out,1);
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

index_y3 = index_y+1;
holder1 = logical(index_y3==(R+1));
index_y3(holder1==1) = 1;
index_x13(holder1==1) = Q-index_x13(holder1==1)+2;
index_x23(holder1==1) = Q-index_x23(holder1==1)+2;
index_x33(holder1==1) = Q-index_x33(holder1==1)+2;

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

function G_out = Quick3DFilter2(im,Gd,G,index)
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

% Set up 2 by 2 matrix per voxel:
A(:,1,1) = a12(:);
A(:,1,2) = a13(:);
A(:,2,1) = a22(:)-e1(:);
A(:,2,2) = a23(:);
% Set up 2 by 1 vector per voxel:
b(:,1) = -a11(:)+e1(:);
b(:,2) = -a12(:);

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
Px_out = reshape(1./r,s1,s2,s3);
Py_out = reshape(coeffs(:,1)./r,s1,s2,s3);
Pz_out = reshape(coeffs(:,2)./r,s1,s2,s3);


% % Set up 2 by 2 matrix per voxel:
% A(:,1,1) = a11(:)-e1(:);
% A(:,1,2) = a12(:);
% A(:,2,1) = a12(:);
% A(:,2,2) = a22(:)-e1(:);
% % Set up 2 by 1 vector per voxel:
% b(:,1) = -a13(:);
% b(:,2) = -a23(:);

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
% Pz_out = reshape(1./r,s1,s2,s3);
% Px_out = reshape(coeffs(:,1)./r,s1,s2,s3);
% Py_out = reshape(coeffs(:,2)./r,s1,s2,s3);


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