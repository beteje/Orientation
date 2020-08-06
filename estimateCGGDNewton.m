%Mike Novey
%University of Maryland Baltimore County
%2/05/09    initial
%%%Function to estimate the shape parameter, c=1 Gaussian, and Covariance
%Of the bivariate GGD distribution using Newton method
%Xa = 2 x N complex augmented matrix Xa = [X, X*]';
%c = 1 is Gausssian
%R is augmented covariance matrix [[a conj(b)]' [b a]']
%Submitted to IEEE Transactions on SIgnal Processing 2009

function [c,Ra, count] = estimateCGGDNewton(Xa)

mom=1;
N = size(Xa,2);

R = cov(Xa'); %Consistant estimator but not MLE
Xr = [real(Xa(1,:))' imag(Xa(1,:))']';
Rr = cov(Xr');
bestC = max(.5,.1*randn+1); %start at Gaussian
%Moment estimator to get c in ball park prior to MLE
if mom == 1
    %Moment estimator for c to get in ball park
    bestVal = 1e10;
    x = real(Xa(1,:));y = imag(Xa(1,:));
    mm = mean(x.^4)/mean(x.^2)^2 + mean(y.^4)/mean(y.^2)^2;
    for cc = [[.1:.1:.45] [.5:.5:1.5] [2:4]]
        temp = 3*gamma(1/cc)*gamma(3/cc)/gamma(2/cc)^2 -mm; %root
        if abs(temp) < bestVal
            bestVal = abs(temp);
            bestC = cc;
        end
    end
end
c = bestC;
%%%%%%%test
% Rr  = [[.6 .22]' [.22 .4]'];
% R  = [[(Rr(1,1)+Rr(2,2)) (Rr(1,1)-Rr(2,2) + j*Rr(1,2)*2)]' [(Rr(1,1)-Rr(2,2) - j*Rr(1,2)*2) (Rr(1,1)+Rr(2,2))]'];
% X = zeros(2,2);
% X(1,1) = .3-j*.22;
% X(2,1) = .3+j*.22;
% X(1,2) = .3-j*.22;
% X(2,2) = .3+j*.22;
% 
% Xr = [real(X(1,:))' imag(X(1,:))']';
% c = .75
% N = size(X,2);
%%%%%%%

paramVec = [Rr(1,1) Rr(2,2) Rr(1,2) c]';


oldParam = zeros(4,1);
maxCount = 25;
count = 0;
delt = 100;

while count < maxCount && delt > 1e-5
    %put paramvec into variables
    Rr(1,1) = paramVec(1);
    Rr(2,2) = paramVec(2);
    Rr(1,2) = paramVec(3);
    Rr(2,1) = paramVec(3);
    c = (paramVec(4));

    %Do complex one also
    aaa = paramVec(1,1) + paramVec(2,1);
    bbb = paramVec(1,1) - paramVec(2,1) + j*paramVec(3,1)*2;
    R = [[aaa (bbb)]' [conj(bbb) aaa]'];

    beta = c*gamma(2/c)/(pi*gamma(1/c)^2);
    eta = gamma(2/c)/(2*gamma(1/c));
    c2p = log(eta) - inv(c)*(2*psi(2/c)-psi(1/c));
    np = psi(1/c)*eta/c^2 - 2*eta*psi(2/c)/c^2;
    %R = [[2 .25]' [.25 3]'];
    %Trace runs slower
    %     xRxC = real(trace((X'*inv(R)*X).^c));
    %     dirXRX = real(trace( log(X'*inv(R)*X).* ((X'*inv(R)*X).^c)));
    %     dirXRX2 = real(trace( (log(X'*inv(R)*X).^2).* ((X'*inv(R)*X).^c)));
    xRxC = 0;
    dirXRX = 0;
    dirXRX2 = 0;
    da= 0;
    daa1 = 0;
    db= 0;
    dbb1 = 0;
    dr = 0;
    drr1 = 0;
    dab1 = 0;
    dar1 = 0;
    dbr1 = 0;
    dac1 = 0;
    dac2 = 0;
    dbc1 = 0;
    dbc2 = 0;
    drc1 = 0;
    drc2 = 0;

    detRr = det(Rr);
    for n = 1:N
        temp = (Xa(:,n)'*inv(R)* Xa(:,n));
        xRxC = xRxC +   real(temp^(c));
        dirXRX = dirXRX + real(log(temp)*temp^(c));
        dirXRX2 = dirXRX2 + real(log(temp)^2*temp^(c));
        xrx = (Xr(:,n)'*inv(Rr)*Xr(:,n));
       
        %%a
        ma = ((Xr(2,n)^2  - Rr(2,2)*xrx))*inv(detRr);
        
        da = da + (xrx)^(c-1)*ma;
        
        maa = (-Xr(2,n)^2 * Rr(2,2) -Rr(2,2)*ma*detRr+Rr(2,2)^2*xrx)/detRr^2;
        temp1 = (c-1)*xrx^(c-2)*ma^2 + xrx^(c-1)*maa;
        daa1 = daa1  + (temp1) ;
        
        %%%b
        mb = ((Xr(1,n)^2  - Rr(1,1)*xrx))*inv(detRr);
        db = db +  (xrx)^(c-1)*mb ;
        
        mbb = (-Xr(1,n)^2 * Rr(1,1) -Rr(1,1)*mb*detRr+Rr(1,1)^2*xrx)/detRr^2;
        temp1 = (c-1)*xrx^(c-2)*mb^2 + xrx^(c-1)*mbb;
        dbb1 = dbb1  + (temp1) ;

        %%% rho
        rho = Rr(1,2);
        mr = (( rho*xrx)- Xr(1,n)*Xr(2,n))*inv(detRr)*2;
        dr = dr + (xrx)^(c-1)*mr;
        
        mrr = 2*inv(detRr^2)*( (xrx + rho*mr)*detRr + 2*rho^2*xrx -2 *rho*Xr(1,n)*Xr(2,n));
        temp1 = (c-1)*xrx^(c-2)*mr^2 + xrx^(c-1)*mrr;
        drr1 = drr1  + temp1;
        
        %a b
        mab = (-Xr(2,n)^2 * Rr(1,1) - detRr*(xrx+Rr(2,2)*mb) + Rr(2,2)*Rr(1,1)*xrx)/detRr^2;
        temp1 = (c-1)*xrx^(c-2)*ma*mb + xrx^(c-1)*mab;
        dab1 = dab1  + (temp1) ;

        %a r
        mar = (2*Xr(2,n)^2 * rho - Rr(2,2)*(mr*detRr+2*rho*xrx))/detRr^2;
        temp1 = (c-1)*xrx^(c-2)*ma*mr + xrx^(c-1)*mar;
        dar1 = dar1  + (temp1) ;

        %b r
        mbr = (2*Xr(1,n)^2 * rho - Rr(1,1)*(mr*detRr+2*rho*xrx))/detRr^2;
        temp1 = (c-1)*xrx^(c-2)*mb*mr + xrx^(c-1)*mbr;
        dbr1 = dbr1  + (temp1) ;

         %a c
        dac1 = dac1  + xrx^(c-1)*ma  ;
        dac2 = dac2  + xrx^(c-1)*log(xrx)*ma ;
        
         %b c
        dbc1 = dbc1  + xrx^(c-1)*mb  ;
        dbc2 = dbc2  + xrx^(c-1)*log(xrx)*mb ;
         %r c
        drc1 = drc1  + xrx^(c-1)*mr  ;
        drc2 = drc2  + xrx^(c-1)*log(xrx)*mr ;

    end
   %Gradient
   
    gc = N*(inv(c) - inv(c^2)*2*psi(2/c)+inv(c^2)*2*psi(1/c))-(eta^c)*(c2p*xRxC + dirXRX);
    
    da = -N*.5 * Rr(2,2)/detRr - (eta^c * c * da);
    db = -N*.5 * Rr(1,1)/detRr - (eta^c * c * db);
    dr = N*Rr(1,2)/detRr - (eta^c * c * dr);
    %Second dir
    daa = N*.5*Rr(2,2)^2*inv(detRr^2) - eta^c*c*(daa1);
    dbb = N*.5*Rr(1,1)^2*inv(detRr^2) - eta^c*c*(dbb1);
    drr = N*(detRr + 2*Rr(1,2)^2)*inv(detRr^2)- eta^c*c*drr1;
    dab = -N*.5*(detRr-Rr(2,2)*Rr(1,1))*inv(detRr^2) - eta^c*c*(dab1);
    dar = -N*Rr(2,2)*rho*inv(detRr^2) - eta^c*c*(dar1);
    dbr = -N*Rr(1,1)*rho*inv(detRr^2) - eta^c*c*(dbr1);
    dac = -((c*eta^c*(log(eta)+(c*np/eta)) + eta^c)*dac1 + eta^c*c*dac2); 
    dbc = -((c*eta^c*(log(eta)+(c*np/eta)) + eta^c)*dbc1 + eta^c*c*dbc2); 
    drc = -((c*eta^c*(log(eta)+(c*np/eta)) + eta^c)*drc1 + eta^c*c*drc2); 
    %%Second dir c
    A = N*((4*psi(2/c)/c^3) + (4*psi(1,2/c)/c^4)-(1/c^2) - (4*psi(1/c)/c^3) - (2*psi(1,1/c)/c^4));
    %Dir c2^c
    dc2C = log(eta)*(eta^c) - c*(eta^(c-1))*(eta*2*psi(2/c)/c^2 - eta*psi(1/c)/c^2);
    dc2p= -((psi(1/c) - 2*psi(2/c))/c^2) - ((psi(1,1/c) - 4*psi(1,2/c))/c^3) - ...
        ((2*psi(2/c)/c^2) - psi(1/c)/c^2);
    B = dc2C*c2p *xRxC + eta^c*(dc2p*xRxC + c2p*dirXRX);
    C = dc2C*dirXRX + eta^c*(dirXRX2);
    ggc = A - B - C; %Testted with MAthacd
    
    %%%%%%%%%%%%%%%%%%%Find other dirs
    H = [[daa dab dar dac]' [dab dbb dbr dbc]' [dar dbr drr drc]' [dac dbc drc ggc]'];
    D = [da db dr gc]';
    
    oldParam = paramVec;
    testVal = rcond(H);
    if isnan(testVal) || abs(testVal) < 1e-7 || isinf(testVal)
        count = maxCount;
        break;
    end
    paramVec = paramVec - inv(H)*D;
    
    paramVec(4) = min(5,max(.1,real(paramVec(4)))); %Newton update with no
    %negatives
    % delt = norm(Rold-R);
    delt = norm(oldParam-paramVec);
     %disp([num2str(count) ' ' num2str(delt)]);
    count=count+1;
end %while
%Put into complex augmented form
a = paramVec(1,1) + paramVec(2,1);
b = paramVec(1,1) - paramVec(2,1) - j*paramVec(3,1)*2;
Ra = [[a conj(b)]' [b a]'];
%R = [[paramVec(1,1) paramVec(3,1)]' [paramVec(3,1) paramVec(2,1)]'];
c = paramVec(4,1);

