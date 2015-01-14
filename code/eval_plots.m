function eval_plots(w1,w2,rho)

Amat = @(k) [-w1 0 1+k.^2 1+k.^2; 0 0 0 1-2.*k.^2; 1 0 w1.*k.^2 w1.*(1+k.^2); 1./rho 1-1./rho -w1.*k.^2./rho w2+(w1/rho+2*w2).*k.^2];

Kvals = (0:.01:100);

Evals = zeros(length(Kvals),4);

for jj=1:length(Kvals)
   
    Evals(jj,:) = eig(Amat(Kvals(jj)));
    
end

clf

figure(1)

hold on

for jj=1:4
   
    plot(Kvals,real(Evals(:,jj)))
    
end

hold off

figure(2)

hold on

for jj=1:4
   
    plot(Kvals,imag(Evals(:,jj)))
    
end

hold off