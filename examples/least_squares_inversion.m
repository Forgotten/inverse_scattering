%% setup scaling parameters
clear
addpath('../src')  
%% butterfly and helmholtz parameters

% number of Gaussian inclusions

nn = 2;
freq = 10; % Hz
omega = 2*pi*freq;

delta_m = 0.2; % amplitude
%     delta_m = 0
lbd = 1/(sqrt(1+delta_m)*freq);

L = 4;
s = 5; % make leaf node several wavelengths?
N = (2^L)*s;
h =  1/(N-1);

fprintf('freq = %.2f Hz\n', freq)
fprintf('lbd  = %.2e \n', lbd)
fprintf('h = %.2e \n', h)

fprintf('N = %d\n', N)
fprintf('PPW = %d\n', floor(lbd/h))

%% setup the model and the domain
% background wavespeed
% N = 100;
c = ones(N,N);

% size of the model in interior domain
nxi  = size(c,2);
nyi  = size(c,1);
ni   = nxi*nyi;

xi = h*(0:nxi-1) - 0.5;
yi = h*(0:nyi-1) - 0.5;

[Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
npml = 20;
nx = nxi + 2*npml;
ny = nyi + 2*npml;
n  = nx*ny;

x  = [xi(1)+(-npml:-1)*h xi xi(end)+(1:npml)*h];
y  = [yi(1)+(-npml:-1)*h yi yi(end)+(1:npml)*h];

[X,Y]   = meshgrid(x,y);

% order of accuracy
order = 8;

% intensity of the pml absorbtion
sigmaMax = 80;

%%  we create a structure (this makes our lives easier)
% this should be encapsulated but for now we just use it.

% putting together all the different properties
properties.nx = nx;
properties.ny = ny;

properties.nxi = nxi;
properties.nyi = nyi;

properties.npml = npml;
properties.omega = omega;

properties.N = N;

properties.X = X;
properties.Y = Y;

% intensity of the pml absorbtion
properties.sigmaMax = sigmaMax;
properties.h = h;

% order of accuracy
% we use a different order to avoid the inverse crime
properties.order = 4;

%%  We generate the data 

% normpdf with amplitude scale = 1 fixed
normal = @(x, mu, sigma)(normpdf(x, mu, sigma)*sqrt(2*pi)*sigma);

% number of scatterers (uniform between [2,4])

eta = (0.2*normpdf(Xi(:),0.5,5*h).*normpdf(Yi(:),0.5,5*h));
centres = 0.4*rand(nn,2) -0.2 ;

sigma_h = 2*h;

for ii = 1:nn
    if ii == 1
        eta = delta_m*normal(Xi(:),centres(ii,1),sigma_h).*...
                      normal(Yi(:),centres(ii,2),sigma_h);
    else
        eta = eta+delta_m*normal(Xi(:),centres(ii,1),sigma_h).*...
                          normal(Yi(:),centres(ii,2),sigma_h);
    end
end
% extend the model to the simulation domain
m = 1 + eta;
eta_ext=ExtendModel(eta,nxi,nyi,npml);
mext=ExtendModel(m,nxi,nyi,npml);

figure(1); clf();
DisplayField(1./sqrt(m),xi,yi);
title('Velocity');

% Helmholtz matrix without auxiliaries and with analytic 
%  differentiation of the pml profile
H1=HelmholtzMatrix(mext,nx,ny,npml,h,...
                   sigmaMax,order,omega,'compact_explicit');


%% computing the farfield data

% computing the number of directions
% Ntheta = 100;
Ntheta = N;
dtheta = 2*pi/(Ntheta);

theta = linspace(pi, 3*pi-dtheta, Ntheta);
d = [cos(theta).' sin(theta).'];

% building the right hand sides
U_in =  exp(1i*omega*(X(:)*d(:,1).'+ Y(:)*d(:,2).'));
S    = bsxfun(@times, -omega^2*eta_ext, U_in);

% solving the equation
tic;
U = H1\S;
t_f = toc;

% printing the time of the solution
fprintf('Time elapsed of the computation = %.4e [s]\n',t_f );


% sampling the wavefield in the interpolated data points
scatter = zeros(Ntheta, Ntheta);
theta_r = linspace(0, 2*pi-dtheta, Ntheta);
r = [cos(theta_r).' sin(theta_r).'];

points_query = 0.5*r;

Projection_mat = zeros(Ntheta, nx, ny);

for ii = 1:nx
    for jj = 1:ny
        mat_dummy = zeros(nx,ny);
        mat_dummy(ii,jj) = 1;
        Projection_mat(:,ii,jj) = interp2(x,y,...
                                          reshape(mat_dummy, nx, ny),...
                                          points_query(:,1),...
                                          points_query(:,2));
    end
end

% we save the projection matrix into the properties structure
properties.project_mat = sparse(reshape(Projection_mat, Ntheta, nx*ny));

% this is our "real data"
scatter = reshape(Projection_mat, Ntheta, nx*ny)*U;

% define the misfit function (it provides both the misfit and the 
% gradient

eta_0 = 0*eta;
% regularization coefficient
eps_reg = 0.001;

[ l2err, Df] = misfit(scatter, eta_0, properties);
fprintf("l2err is %.4e \n", l2err)
JJ  = @(x) J_starJ(eta_0, x, properties, eps_reg);
tol = 1e-3;
tic;
delta_eta =  gmres(JJ,-Df,10,tol,10);
toc

figure(5);
clf();
subplot(1,3,1)
imagesc(reshape(eta, nxi, nyi))
colorbar();
title('Exact')
subplot(1,3,2)
imagesc(reshape(delta_eta, nxi, nyi))
colorbar();
title('Reconstruction')
subplot(1,3,3)
imagesc(reshape(delta_eta-eta, nxi, nyi))
colorbar();
title('Error')

% % 
% for i = 1:4
% 
%     [ l2err, Df] = misfit(scatter, eta_0, properties);
%     fprintf("l2err is %.4e \n", l2err)
%     JJ  = @(x) J_starJ(eta_0, x, properties);
%     tol = 1e-3;
%     tic;
%     delta_eta =  gmres(JJ,-Df,4,tol,10);
%     toc
%     eta_0 = eta_0 + delta_eta;
% 
% end 
% 
% figure(6);
% clf();
% subplot(1,3,1)
% imagesc(reshape(eta, nxi, nyi))
% colorbar();
% title('Exact')
% subplot(1,3,2)
% imagesc(reshape(eta_0, nxi, nyi))
% colorbar();
% title('Reconstruction')
% subplot(1,3,3)
% imagesc(reshape(eta_0-eta, nxi, nyi))
% colorbar();
% title('Error')
