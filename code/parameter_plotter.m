function parameter_plotter(rho,B)

    omega = (-2:.005:2);
    
    anl = @(lambda1,w1,w2) 2.*(-(1/2)*w1*rho^2.*lambda1.^5 + ( w1*rho^2*w2 + rho^2 + rho./2).*lambda1.^4 + (w1.*rho./2 - 2.*rho.^2.*w2 + rho.^2.*w1 - rho.^2.*w1.*w2.^2./2).*lambda1.^3+...
                            (-3.*rho.^2./2 + rho.^2.*w2.^2 - rho.^2.*w1.*w2 + 3.*rho./2).*lambda1.^2 + (-3.*rho.*w2./2 + rho.^2.*w2 - rho.^2.*w1./2 + w1.*rho).*lambda1...
                            + 1-rho + rho.^2.*w2.^2./2 - rho.*w1.*w2./2 + w2.*(-rho+1).*(2.*rho+1)./(2.*lambda1) + (1/2)*(-rho+1)^2./lambda1.^2);
    
    ad = @(lambda1,w1,w2) ((rho*w1-2*rho*w2-w1).*lambda1.^4+(rho*w1^2-3*rho*w1*w2+6*B+rho-3).*lambda1.^3+...
                          (-rho*w1^2*w2+w1^3+6*B*w1- rho*w1+3*rho*w2-2*w1).*lambda1.^2+(-rho*w1^2+2*rho*w1*w2-w1^2-6*B- rho+1).*lambda1+2*rho*w1-2*w1)./(6.*lambda1);
    
    at = @(lambda1,w1,w2) 2.*lambda1.^3+lambda1.^2.*rho.*(w1-w2)+rho.*(w1-w2)+2.*(1-rho)./lambda1;
    
    avals = zeros(length(omega),length(omega),4,3);
    
    balance = zeros(length(omega),length(omega),4);
    
    evals = zeros(length(omega),length(omega),4);
    
    for jj=1:length(omega)

        for kk=1:length(omega)
            
            w1 = omega(kk); w2 = omega(jj);
            
            rloc = roots([1 -(w2-w1) -(2+w1*w2) (w2-w1) 1-1/rho]);
            rloc = sort(rloc);
            
            for ll=1:4
               
                avals(jj,kk,ll,:) = [anl(rloc(ll),w1,w2); ad(rloc(ll),w1,w2); at(rloc(ll),w1,w2)];
                
                balance(jj,kk,ll) = -log10(abs(anl(rloc(ll),w1,w2)./(ad(rloc(ll),w1,w2).*at(rloc(ll),w1,w2))));
                
                evals(jj,kk,ll) = rloc(ll);
                
            end
               
        end
        
    end
    
% Plot the balances between the different coefficients    


    figure(1)

    subplot(2,2,1)
    
    surf(omega,omega,log10(abs(avals(:,:,1,1))),'LineStyle','none'); view(2);
    
    subplot(2,2,2)
    
    surf(omega,omega,log10(abs(avals(:,:,2,1))),'LineStyle','none'); view(2);
    
    subplot(2,2,3)
    
    surf(omega,omega,log10(abs(avals(:,:,3,1))),'LineStyle','none'); view(2);
    
    subplot(2,2,4)
    
    surf(omega,omega,log10(abs(avals(:,:,4,1))),'LineStyle','none'); view(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(2)

    subplot(2,2,1)
    
    surf(omega,omega,log10(abs(avals(:,:,1,2))),'LineStyle','none'); view(2);
    
    subplot(2,2,2)
    
    surf(omega,omega,log10(abs(avals(:,:,2,2))),'LineStyle','none'); view(2);
    
    subplot(2,2,3)
    
    surf(omega,omega,log10(abs(avals(:,:,3,2))),'LineStyle','none'); view(2);
    
    subplot(2,2,4)
    
    surf(omega,omega,log10(abs(avals(:,:,4,2))),'LineStyle','none'); view(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(3)

    subplot(2,2,1)
    
    surf(omega,omega,log10(abs(avals(:,:,1,3))),'LineStyle','none'); view(2);
    
    subplot(2,2,2)
    
    surf(omega,omega,log10(abs(avals(:,:,2,3))),'LineStyle','none'); view(2);
    
    subplot(2,2,3)
    
    surf(omega,omega,log10(abs(avals(:,:,3,3))),'LineStyle','none'); view(2);
    
    subplot(2,2,4)
    
    surf(omega,omega,log10(abs(avals(:,:,4,3))),'LineStyle','none'); view(2);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(4)
    
    subplot(2,2,1)
    
    surf(omega,omega,balance(:,:,1),'LineStyle','none'); view(2);

    subplot(2,2,2)
    
    surf(omega,omega,balance(:,:,2),'LineStyle','none'); view(2);

    subplot(2,2,3)
    
    surf(omega,omega,balance(:,:,3),'LineStyle','none'); view(2);

    subplot(2,2,4)
    
    surf(omega,omega,balance(:,:,4),'LineStyle','none'); view(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(5)
    
    subplot(2,2,1)
    
    contourf(omega,omega,evals(:,:,1),10,'ShowText','on')

    subplot(2,2,2)
    
    contourf(omega,omega,evals(:,:,2),10,'ShowText','on')

    subplot(2,2,3)
    
    contourf(omega,omega,evals(:,:,3),10,'ShowText','on')

    subplot(2,2,4)
    
    contourf(omega,omega,evals(:,:,4),10,'ShowText','on')

    
% Generate contour plots of all possible coefficients    
%{

    figure(3)
    
    subplot(4,3,1)
    
    contourf(omega,omega,avals(:,:,1,1),5,'ShowText','on')
   
    subplot(4,3,2)
    
    contourf(omega,omega,avals(:,:,1,2),5,'ShowText','on')
    
    subplot(4,3,3)
    
    contourf(omega,omega,avals(:,:,1,3),5,'ShowText','on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(4,3,4)
    
    contourf(omega,omega,avals(:,:,2,1),5,'ShowText','on')
    
    subplot(4,3,5)
    
    contourf(omega,omega,avals(:,:,2,2),5,'ShowText','on')
    
    subplot(4,3,6)
    
    contourf(omega,omega,avals(:,:,2,3),5,'ShowText','on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(4,3,7)
    
    contourf(omega,omega,avals(:,:,3,1),5,'ShowText','on')
    
    subplot(4,3,8)
    
    contourf(omega,omega,avals(:,:,3,2),5,'ShowText','on')
    
    subplot(4,3,9)
    
    contourf(omega,omega,avals(:,:,3,3),5,'ShowText','on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(4,3,10)
    
    contourf(omega,omega,avals(:,:,4,1),5,'ShowText','on')
    
    subplot(4,3,11)
    
    contourf(omega,omega,avals(:,:,4,2),5,'ShowText','on')
    
    subplot(4,3,12)
    
    contourf(omega,omega,avals(:,:,4,3),5,'ShowText','on')   
%}  
end