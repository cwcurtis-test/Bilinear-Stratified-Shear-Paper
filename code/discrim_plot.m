function discrim_plot(rho)

strat = 1 - 1/rho;

omega = (-2:.05:2);

disc = @(w1,w2) 256*strat^3+(-27*w1.^4-36*w1.^3.*w2-2*w1.^2.*w2.^2-36.*w1.*w2.^3-27.*w2.^4-96.*w1.^2-320.*w1.*w2-96.*w2.^2-512).*strat^2 +...
                (4.*w1.^5.*w2.^3+8.*w1.^4.*w2.^4+4.*w1.^3.*w2.^5+18.*w1.^5.*w2+32.*w1.^4.*w2.^2+28.*w1.^3.*w2.^3+32.*w1.^2.*w2.^4+18.*w1.*w2.^5+30.*w1.^4+104.*w1.^3.*w2+116.*w1.^2.*w2.^2+104.*w1.*w2.^3+30.*w2.^4+64.*w1.^2+384.*w1.*w2+64.*w2.^2+256)*strat + ...
                 w1.^6.*w2.^2-2.*w1.^4.*w2.^4+w1.^2.*w2.^6+4.*w1.^6-2.*w1.^5.*w2-4.*w1.^4.*w2.^2+4.*w1.^3.*w2.^3-4.*w1.^2.*w2.^4-2.*w1.*w2.^5+4.*w2.^6+13.*w1.^4-4.*w1.^3.*w2-18.*w1.^2.*w2.^2-4.*w1.*w2.^3+13.*w2.^4+32.*w1.^2-64.*w1.*w2+32.*w2.^2;
             
dvals = zeros(length(omega),length(omega));
    
for jj=1:length(omega)
        
    dvals(jj,:) =  disc(omega(jj),omega);
               
end

disp(min(min(dvals)))

figure(1)

contourf(omega,omega,dvals,10,'ShowText','on')
