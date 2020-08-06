%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference
%Mike Novey and T. Adali, "Circularity and Gaussianity detection
%     using the complex generalized Gaussian distribution" in
% (To appear ) IEEE Signal Processing Letters.,
% June 21, 2009
% Function to detect noncircularity and non-Gaussianity 
% where X1 is input vector and pfa is prob. false alarm
% Returns the detection results for noncircularity (nonCirc),
%   non-Gaussianity (nonG), and non white Gaussian (nonWgn) as well
% as the shape parameter estimate cets
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%cest = 1 Gaussian
function [nonCirc, nonG,nonWgn,cest] = circularityDetector(X1, pfa)
 
%calculate thresholds
threshXchi1 = chi2inv(1-pfa,1);
threshXchi2 = chi2inv(1-pfa,2);
threshXchi3 = chi2inv(1-pfa,3);
N = size(X1,1);
%normalize
%X2 = X1/std(X1);
X = [X1 conj(X1)]';
 
%[cest,C1] = estimateGGDCovShape(X);
[cest,C1] = estimateCGGDNewton(X);
c2 = (gamma(2/cest)/(2*gamma(1/cest)));
c1 = (cest*gamma(2/cest))/(gamma(1/cest)^2);
 
C0UC = [[C1(1,1) 0]' [0 C1(1,1)]']; %Tc only
 
vv1 = (X(1,:)*X(1,:)')/length(X);
C0C = [[vv1 0]' [0 vv1]']; %Twg
 
CS = (X*X')/length(X);
 
TC = 0;
TWG = 0;
TG = 0;
for nn = 1:N
    useC0 = (X(:,nn)'*inv(C0UC)*X(:,nn))^cest;
 
    useC1 = (X(:,nn)'*inv(C1)*X(:,nn))^cest;
    useC0WG = (X(:,nn)'*inv(C0C)*X(:,nn));
    TC = TC + 2*(c2^cest)*(useC0 - useC1);
    TWG = TWG + useC0WG - 2*(c2^cest)*useC1;
    TG = TG  - 2*(c2^cest)*useC1;
end
TC = real(TC -N*log(det(C1)) + N*log(det(C0UC)));
 
TWG = real(TWG + 2*N*log(c1)- N*log(det(C1)) +N*log(det( C0C)));
 
TG = real(TG + 2*N*log(c1) - N*log(det(C1)) +N*log(det(CS)) +2*N);
 
%%%%%%%%bayesian
% uStart = .2;  %remember 1 is Gaussian here
% uStop = 2;
% NumB = 30;
% pB = linspace(uStart,uStop,NumB)';
% dp = pB(2,1) - pB(1,1);
% pU = 1/(uStop-uStart);
% %pU = inv(sqrt(2*pi*.5))*exp(-inv(2*.5)*(pB-2).^2);
% c1B = (2*pB.*gamma(2./pB))./(gamma(1./pB).^2);
% c2B =  (gamma(2./pB)./(2*gamma(1./pB))).^(pB);
% expB = zeros(NumB,1);;
% for bb = 1:NumB
%     expB(bb,1) = real(exp(-c2B(bb,1).*(trace((X'*inv(R1)*X).^pB(bb,1)))));
% end
% numErator = sum(pU.*(c1B/sqrt(abs(det(R1)))).^N .* expB*dp);
% %log(numErator)
% expB = zeros(NumB,1);
% for bb = 1:NumB
%     expB(bb,1) = real(exp(-c2B(bb,1).*(trace((X'*X).^pB(bb,1)))));
% end
% denOm =  sum(pU.*(c1B).^N .* expB)*dp;
% %log(denOm)
% TB = abs(2*(log(real(numErator))-log(real(denOm))));
 
nonCirc = 1.0*(TC > threshXchi2);
nonWgn = 1.0*(TWG > threshXchi3);
nonG = 1.0*(TG > threshXchi1);
