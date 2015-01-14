function movie_plot = kdv_solver_etdrk4(K,Llx,tf,dt,ep,cnl,cd,ct,u0,inter,lam)

% This is a means for solving KdV using the ETDRK4 time integrator.  

KT = 2*K;

dx = Llx/K;

tstep = round(ep*tf/dt);

Dx = 1i.*pi/Llx*[0:K-1 0 -K+1:-1]';
Dop = -Dx.^3*cd/ct;

disp(max(abs(dt*Dop)))

% Build ETDRK4 machinery

    Em = exp(dt*Dop);
    Em2 = exp(dt*Dop/2);
    Mstep = 32;
    r = exp(1i*pi*((1:Mstep)-.5)/Mstep);
    LR = dt*Dop(:,ones(Mstep,1))+r(ones(KT,1),:);

    f1 = dt*real(mean( (exp(LR/2)-1)./LR  ,2)); 
    f2 = dt*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3  ,2)); 
    f3 = dt*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3  ,2)); 
    Q = dt*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3  ,2));

% Transform initial condition    
    
    u = fft(u0); 

no_of_plots = round(tstep/inter); % Number of plots in movie
movie_plot = zeros(no_of_plots,KT); % Movie storage

plot_count = 1; % Current frame 

% Solve KdV equation in time. 

for jj=1:tstep
    
    tc = (jj-1)*dt;
      
    Nu = -Dx.*fft( real(ifft(u)).^2 )*cnl/(2*ct);
    
    a = Em2.*u + Q.*Nu;
    
    Na = -Dx.*fft( real(ifft(a)).^2 )*cnl/(2*ct);
    
    b = Em2.*u + Q.*Na;
    
    Nb = -Dx.*fft( real(ifft(b)).^2 )*cnl/(2*ct);
    
    c = Em2.*a + Q.*(2*Nb-Nu);
    
    Nc = -Dx.*fft( real(ifft(c)).^2 )*cnl/(2*ct);
    
    u = Em.*u + f1.*Nu +2*f2.*(Na+Nb) + f3.*Nc;
    
    if( mod(jj,inter) == 0)  % Capture perturbation at every inter number of time steps
        
            uc = real(ifft(u));
        
            shift = floor(lam*tc/(dx*ep));
            
            uc = circshift(uc,shift);
            
            movie_plot(plot_count,:) = uc; % Track perturbation at time tcurr
            plot_count = plot_count + 1;
        
    end            
    
end

end



