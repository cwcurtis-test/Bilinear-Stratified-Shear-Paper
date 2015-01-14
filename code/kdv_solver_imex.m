function movie_plot = kdv_solver_imex(K,Llx,tf,dt,ep,cnl,cd,ct,u0,inter,lam)

% This is a means for solving the KdV equation using an IMEX scheme that is 
% second order implicit in the linear part and fourth order explicit in the linear part

KT = 2*K;

dx = Llx/K;

tstep = round(ep*tf/dt);

% Build ETDRK4 machinery to generate initial time steps for IMEX scheme.

Dx = 1i.*pi/Llx*[0:K-1 0 -K+1:-1]';
Dop = -Dx.^3*cd/ct;

Lop = (ones(KT,1)-3*dt*Dop/4).^(-1);

Em = exp(dt*Dop);
Em2 = exp(dt*Dop/2);
Mstep = 32;
r = exp(1i*pi*((1:Mstep)-.5)/Mstep);
LR = dt*Dop(:,ones(Mstep,1))+r(ones(KT,1),:);

f1 = dt*real(mean( (exp(LR/2)-1)./LR  ,2)); 
f2 = dt*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3  ,2)); 
f3 = dt*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3  ,2)); 
Q = dt*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3  ,2));

no_of_plots = round(tstep/inter); % Number of plots in movie
movie_plot = zeros(no_of_plots,KT); % Movie storage

plot_count = 1; % Current frame 

u = fft(u0); 
uprev = zeros(KT,4);
uprev(:,1) = u;

% Run ETDRK4 in order to generate enough information to start IMEX scheme.  

for jj=1:3
    
    Nu = -Dx.*fft( real(ifft(u)).^2 )*cnl/(2*ct);
    
    a = Em2.*u + Q.*Nu;
    
    Na = -Dx.*fft( real(ifft(a)).^2 )*cnl/(2*ct);
    
    b = Em2.*u + Q.*Na;
    
    Nb = -Dx.*fft( real(ifft(b)).^2 )*cnl/(2*ct);
    
    c = Em2.*a + Q.*(2*Nb-Nu);
    
    Nc = -Dx.*fft( real(ifft(c)).^2 )*cnl/(2*ct);
    
    u = Em.*u + f1.*Nu +2*f2.*(Na+Nb) + f3.*Nc;
    
    uprev(:,jj+1) = u;
 
end

u1 = uprev(:,1);
u2 = uprev(:,2); 
u3 = uprev(:,3); 
u4 = uprev(:,4); 

% Run IMEX scheme 

for jj=5:tstep
    
    
    tc = (jj-1)*dt;
    
    u5 = Lop.*(u4 + dt*Dop.*u3/4 - dt*Dx.*cnl.*fft(-9*real(ifft(u1)).^2 + 37*real(ifft(u2)).^2 - 59*real(ifft(u3)).^2 + 55*real(ifft(u4)).^2)./(48.*ct));
    
    u1 = u2;
    u2 = u3;
    u3 = u4;
    u4 = u5;
    
    if( mod(jj,inter) == 0)  % Capture perturbation at every inter number of time steps
        
            uc = real(ifft(u5));
            
            shift = floor(lam*tc/(dx*ep));
            
            uc = circshift(uc,shift);
            
            movie_plot(plot_count,:) = uc; % Track perturbation at time tcurr
            plot_count = plot_count + 1;
        
    end            
    
end

end



