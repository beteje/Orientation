function [Circularity_XY,Circularity_YZ,Circularity_XZ,R_length,az_mean,el_mean] = Directionality_measures(az_est,el_est)

% convert of radians:
az_est = az_est*pi./180;
el_est = el_est*pi./180;

% Generate 3D vector:
Vx = sin(el_est).*cos(2*az_est);
Vy = sin(el_est).*sin(2*az_est);
Vz = cos(el_est);

% Generate complex numbers:
alpha = atan2(Vy,Vx);
Vxy = sqrt(Vx.^2+Vy.^2).*exp(1j.*2.*alpha);
alpha = atan2(Vz,Vy);
Vyz = sqrt(Vz.^2+Vy.^2).*exp(1j.*2.*alpha);
Vxz = Vx + 1j.*Vz;

Vxy = Vx + 1j.*Vy;
Vyz = Vy + 1j.*Vz;
Vxz = Vx + 1j.*Vz;

% Calculate R values:
Rxy = mean(abs(Vxy(:)).^2);
Rxz = mean(abs(Vyz(:)).^2);
Ryz = mean(abs(Vxz(:)).^2);

% Calculate R1 values:
R1xy = mean(Vxy(:).^2);
R1xz = mean(Vyz(:).^2);
R1yz = mean(Vxz(:).^2);

% Determine the three circularity numbers:
Circularity_XY = abs(R1xy).^2./Rxy.^2;
Circularity_XZ = abs(R1xz).^2./Rxz.^2;
Circularity_YZ = abs(R1yz).^2./Ryz.^2;


% figure; subplot(1,3,1); scatter(real(Vxy(:)),imag(Vxy(:)),'.'); xlim([-1,1]); ylim([-1,1]); axis square; title(['Circ XY = ',num2str(Circularity_XY,3)]);
% subplot(1,3,2); scatter(real(Vyz(:)),imag(Vyz(:)),'.'); xlim([-1,1]); ylim([-1,1]); axis square; title(['Circ YZ = ',num2str(Circularity_YZ,3)]);
% subplot(1,3,3); scatter(real(Vxz(:)),imag(Vxz(:)),'.'); xlim([-1,1]); ylim([-1,1]); axis square; title(['Circ XZ = ',num2str(Circularity_XZ,3)]);
% 
% figure; subplot(2,1,1); histogram(az_est.*180./pi,linspace(0,180,180/5)); title('Azimuth');
% subplot(2,1,2); histogram(el_est.*180./pi,linspace(0,180,180/5)); title('Elevation');

% Adjust for axial data:
Vx = sin(el_est).*cos(2*az_est);
Vy = sin(el_est).*sin(2*az_est);
Vz = cos(el_est);

% Spherical mean:
x(1) = mean(Vx(:));
x(2) = mean(Vy(:));
x(3) = mean(Vz(:));

% Mean resultant length:
R_length = sqrt(sum(abs(x).^2));

% Mean direction:
x_mean = x./R_length;
az_mean = atan2(x_mean(2),x_mean(1));        
az_mean(az_mean<0) = az_mean(az_mean<0)+2*pi;
az_mean = az_mean.*180./pi/2;
el_mean = acos(x_mean(3))*180/pi;



