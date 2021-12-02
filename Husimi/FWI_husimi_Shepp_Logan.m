%% setup scaling parameters
clear
addpath('../src')  
%% helmholtz parameters

% this script performs a monochromatic inversion of a perturbation with 
% multiscale features. 

% The main objective of this code if to showcase the cycle-skipping issue
% in which the optimization get stuck on a local minimum

sigma = 0.01;

nn = 1;
freq = 10; % Hz
omega = 2*pi*freq;

delta_m = 0.5; % amplitude
%     delta_m = 0
lbd = 1/(sqrt(1+delta_m)*freq);

a = 2;  % size of domain has to be an integer

L = 4;
s = 5; % make leaf node several wavelengths?  % what does this mean...?
N = a*(2^L)*s;
h =  a/(N-1);

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

xi = h*(0:nxi-1) - 0.5*a;
yi = h*(0:nyi-1) - 0.5*a;

[Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
npml = 20;
nx = nxi + 2*npml;
ny = nyi + 2*npml;
n  = nx*ny;

x  = [xi(1)+(-npml:-1)*h,xi,xi(end)+(1:npml)*h];
y  = [yi(1)+(-npml:-1)*h,yi,yi(end)+(1:npml)*h];

[X,Y] = meshgrid(x,y);

% order of accuracy
order = 8;

% intensity of the pml absorbtion
sigmaMax = 80;

% Number of sources/receivers
Ntheta_s = 40;
Ntheta_r = 40;
% Number of directions
Ntheta_i = 20;
Ntheta_o = 20;

%% sources and measurement boundary
xb = 0; yb = 0; rb = 0.4;
% Width
sigma0 = 2^(-2);
Cd = sqrt(2)*(sigma0/sqrt(pi))^(3/2);
% Prototype function
Sfunc = @(x,y,thetax,thetay) omega^(5/2)*Cd*...
            exp(-omega^2*sigma0^2/2*(x.^2+y.^2) + 1i*omega*(thetax.*x+thetay.*y));

%% Gabor filter
phi = @(x,y,thetax,thetay) 0.5*(omega/pi)^(3/2)*...
                exp(-omega/2*(x.^2+y.^2) + 1i*omega*(thetax.*x+thetay.*y));

%%  we create a structure (this makes our lives easier)
% this should be encapsulated later

% putting together all the different properties
properties.nx = nx;
properties.ny = ny;
properties.n = n;

properties.nxi = nxi;
properties.nyi = nyi;

properties.npml = npml;
properties.omega = omega;

properties.Ntheta_s = Ntheta_s;
properties.Ntheta_i = Ntheta_i;
properties.Ntheta_r = Ntheta_r;
properties.Ntheta_o = Ntheta_o;
% properties.Ntheta_t = Ntheta_p*Ntheta_d;

% intensity of the pml absorbtion
properties.sigmaMax = sigmaMax;
properties.h = h;

% order of accuracy
% we use a different order to avoid the inverse crime
properties.order = 4;

%%  We generate the data 

%%%
%%% Generate the medium %%%
%%%

% genaration of a random reference perturbation

% we generate the Shepp-Logan phantom rescaled
P = delta_m*phantom('Modified Shepp-Logan',N/2);

eta = zeros(N, N);

eta(N/4+1:3*N/4,N/4+1:3*N/4) = P;



%% extend the model to the simulation domain

% we define the full slowness squared, the background 
% plus the perturbartion
m = 1 + eta;

% we extend it to conside the PMLs
mext = ExtendModel(m,nxi,nyi,npml);

figure(1); clf();
DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
title('Velocity');


%%%
%%% Generate the Helmholtz matrix %%%
%%%

% Helmholtz matrix without auxiliaries and with analytic 
%  differentiation of the pml profile
H1=HelmholtzMatrix(mext,nx,ny,npml,h,...
                   sigmaMax,order,omega,'compact_explicit');
H1 = H1'; % The Helmholtz code produces adjoint matrix (TO FIX)

                          
%%%
%%% computing scattered field %%%
%%%

% computing the number of directions


dtheta_s = 2*pi/Ntheta_s;
dtheta_i = pi/(Ntheta_i+1);

theta_s = repmat(linspace(pi, 3*pi-dtheta_s, Ntheta_s),...
                 Ntheta_i, 1);

theta_i = repmat(linspace(-pi/2+dtheta_i, pi/2-dtheta_i, Ntheta_i)',...
                 1, Ntheta_s);

theta_i = theta_i + theta_s + pi;   

% computing the directions
dir = [cos(theta_i(:)),... 
       sin(theta_i(:))];

% and the positions of the centers of the beams
pos = [xb + rb*cos(theta_s(:)), ...
       yb + rb*sin(theta_s(:))];

% building the right hand sides
S = Sfunc(X(:)-pos(:,1).', Y(:)-pos(:,2).',...
         dir(:,1).'      , dir(:,2).' );


S_plot = reshape(S(:,1),ny,nx); 
S_plot = S_plot(npml+1:npml+nxi, npml+1:npml+nyi);
figure(99); clf();
DisplayField(S_plot,xi,yi); set(gca,'YDir','normal');
title('Source');


% solving the equation
tic;
U = H1\S;
t_f = toc;

% printing the time of the solution
fprintf('Time elapsed of the computation = %.4e [s]\n',t_f );

U_plot = reshape(U(:,1),ny,nx); 
U_plot = U_plot(npml+1:npml+nxi, npml+1:npml+nyi);
figure(98); clf();
DisplayField(U_plot,xi,yi); set(gca,'YDir','normal');
title('Solution');

% we save the source matrix into the properties structure
properties.S = S;

%%%
%%% Perform Husimi transform %%%
%%%

% sampling the wavefield
dtheta_r = 2*pi/Ntheta_r;
dtheta_o = pi/(Ntheta_o+1);

theta_r = repmat(linspace(0, 2*pi-dtheta_r, Ntheta_r),...
                 Ntheta_o, 1);
theta_o = repmat(linspace(-pi/2+dtheta_o, pi/2-dtheta_o, Ntheta_o)',...
                 1, Ntheta_r);

theta_o = theta_o + theta_r;


dir_r = [cos(theta_o(:)), sin(theta_o(:))];
pos_r = [xb + rb*cos(theta_r(:)),...
         yb + rb*sin(theta_r(:))];

husimi_mat = phi(pos_r(:,1)-X(:).', pos_r(:,2)-Y(:).',...
                 dir_r(:,1)       , dir_r(:,2) );

% we save the husimi matrix into the properties structure
properties.husimi_mat = husimi_mat;


%%%
%%% this is our "real data" %%%
%%%

scatter = abs(h^2*husimi_mat*U).^2;

%% Nonlinear least square
% define the misfit function (it provides both the misfit and the 
% gradient

J = @(x) misfit_husimi(scatter, x, properties);

% J = @(x) misfit_husimi_gpu(scatter, x, properties);

% prepare the options for the optimization loop 
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations', 5,...
    'OptimalityTolerance', 1e-5, ...
    'Display', 'iter-detailed');

% option to check the gradients using FD (it takes very loooong)

% options = optimoptions('fminunc','Algorithm','quasi-newton',...
%                         'SpecifyObjectiveGradient',true,...
%                         'CheckGradients', true,...
%                         'FiniteDifferenceType', 'central',...
%                         'MaxIterations', 10000,...
%                         'Display', 'iter-detailed');
%

% running the optimization routine
[result,fval,exitflag,output] = fminunc(J,0*eta,options);

% plotting the result
figure(5);
clf();
subplot(1,3,1)
imagesc(xi,yi,reshape(eta, nxi, nyi)); pbaspect([1 1 1]); 
set(gca,'YDir','normal'); colorbar();
title('Exact')
subplot(1,3,2)
imagesc(xi,yi,reshape(result, nxi, nyi)); pbaspect([1 1 1]);
set(gca,'YDir','normal'); colorbar();
title('Reconstruction')
subplot(1,3,3)
imagesc(xi,yi,reshape(result-eta, nxi, nyi)); pbaspect([1 1 1]);
set(gca,'YDir','normal'); colorbar();
title('Error')

set(gcf, 'Position',  [100, 100, 1500, 400])


saveas(gcf,fullfile('..','plots',...
    ['n_reconstruct_shepp_logan_','_An',num2str(delta_m),...
        '_k',num2str(omega),'_h',num2str(h),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai',int2str(Ntheta_i),'_Nthetas',int2str(Ntheta_s),...
        '_Nthetao',int2str(Ntheta_o),'_Nthetar',int2str(Ntheta_r),'.eps']),'epsc')

% % print and save
% print(fullfile('..','plots',...
%     ['n_reconstruct_medium_delocalized_',int2str(sigma),'_An',num2str(delta_m),...
%         '_k',num2str(omega),'_h',num2str(h),...%'_a',num2str(a),'_b',num2str(b),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai',int2str(Ntheta_i),'_Nthetas',int2str(Ntheta_s),...
%         '_Nthetao',int2str(Ntheta_o),'_Nthetar',int2str(Ntheta_r),'.pdf']),'-dpdf');

savefig(fullfile('..','plots',...
    ['n_reconstruct_shepp_logan_','_An',num2str(delta_m),...
        '_k',num2str(omega),'_h',num2str(h),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai',int2str(Ntheta_i),'_Nthetas',int2str(Ntheta_s),...
        '_Nthetao',int2str(Ntheta_o),'_Nthetar',int2str(Ntheta_r),'.fig']));
              
save(fullfile('..','data',...
    ['n_reconstruct_shepp_logan_','_An',num2str(delta_m),...
        '_k',num2str(omega),'_h',num2str(h),...%'_a',num2str(a),'_b',num2str(b),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai',int2str(Ntheta_i),'_Nthetas',int2str(Ntheta_s),...
        '_Nthetao',int2str(Ntheta_o),'_Nthetar',int2str(Ntheta_r),'.mat']),...
        'eta','result');

