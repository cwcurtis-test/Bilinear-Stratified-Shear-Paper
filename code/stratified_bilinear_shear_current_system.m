function stratified_bilinear_shear_current_system(K,Llx,tf,w1,w2,rho,B,ep)

% This program solves all four different KdV equations and plots the
% computed surface and internal wave heights.  

KT = 2*K;

dx = Llx/K;

dt = 5e-4; tstep = round(ep*tf/dt);

X = (-Llx:dx:Llx-dx)';

rloc = sort(roots([1 -(w2-w1) -(2+w1*w2) (w2-w1) 1-1/rho]));

display('Wave numbers')
display(rloc)

% Build functions to compute various coefficients to use in KdV equations 

anl = @(lambda1,w1,w2) 2.*(-(1/2)*w1*rho^2.*lambda1.^5 + ( w1*rho^2*w2 + rho^2 + rho./2).*lambda1.^4 + (w1.*rho./2 - 2.*rho.^2.*w2 + rho.^2.*w1 - rho.^2.*w1.*w2.^2./2).*lambda1.^3+...
                            (-3.*rho.^2./2 + rho.^2.*w2.^2 - rho.^2.*w1.*w2 + 3.*rho./2).*lambda1.^2 + (-3.*rho.*w2./2 + rho.^2.*w2 - rho.^2.*w1./2 + w1.*rho).*lambda1...
                            + 1-rho + rho.^2.*w2.^2./2 - rho.*w1.*w2./2 + w2.*(-rho+1).*(2.*rho+1)./(2.*lambda1) + (1/2)*(-rho+1)^2./lambda1.^2);
    
ad = @(lambda1,w1,w2) ((rho*w1-2*rho*w2-w1).*lambda1.^4+(rho*w1^2-3*rho*w1*w2+6*B+rho-3).*lambda1.^3+...
                          (-rho*w1^2*w2+w1^3+6*B*w1- rho*w1+3*rho*w2-2*w1).*lambda1.^2+(-rho*w1^2+2*rho*w1*w2-w1^2-6*B- rho+1).*lambda1+2*rho*w1-2*w1)./(6.*lambda1);
    
at = @(lambda1,w1,w2) (2.*lambda1.^3 + lambda1.^2.*rho.*(w1-w2) + rho.*(w1-w2) + 2.*(1-rho)./lambda1);

% Pick initial surface and internal wave heights.  Likewise choose surface and internal velocity potentials. 

n1 = .5.*sech(X).^2;
n2 = .5.*sech(X).^2; 
q1 = .5.*sech(X).^2; 
q2 = .5.*sech(X).^2;

% Compute initial conditions to use in KdV equations

u0 = zeros(KT,4);

for kk=1:KT
    u0(kk,:) = ((at(rloc,w1,w2)).^(-1).*(rloc.*n1(kk)+((w1+rloc).*rloc-1).*((rho-1).*n2(kk)./rloc + rho.*q2(kk)) + q1(kk)))';
end

% Compute and display the amplitudes of the initial conditions that go into the parameter independent KdV equations.  

bal = anl(rloc,w1,w2)./(6.*ad(rloc,w1,w2));

display('Amplitudes of initial conditions for rescaled KdV')
disp([min(bal(1)*u0(:,1)),max(bal(1)*u0(:,1))])
disp([min(bal(2)*u0(:,2)),max(bal(2)*u0(:,2))])
disp([min(bal(3)*u0(:,3)),max(bal(3)*u0(:,3))])
disp([min(bal(4)*u0(:,4)),max(bal(4)*u0(:,4))])

inter = 15; % Number of time steps between frames of movie
no_of_plots = round(tstep/inter); % Number of plots in movie
movie_plot = zeros(no_of_plots,KT,2); % Movie storage

movie_storage = zeros(no_of_plots,KT,4);

% Solve each KdV equation 

for jj=1:4
   
    movie_storage(:,:,jj) = kdv_solver_imex(K,Llx,tf,dt,ep,anl(rloc(jj),w1,w2),ad(rloc(jj),w1,w2),at(rloc(jj),w1,w2),u0(:,jj),inter,rloc(jj));
        
end
    
% Compute surface and internal wave displacements.  

for jj=1:4

    movie_plot(:,:,1) = movie_plot(:,:,1) + (1-rho-rho*w2*rloc(jj)+rho*rloc(jj)^2)*movie_storage(:,:,jj);
    movie_plot(:,:,2) = movie_plot(:,:,2) + movie_storage(:,:,jj);

end
    
Movie_Maker_1d(movie_plot,X,no_of_plots,'two_vortex_patches')

end

