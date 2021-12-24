%% setup scaling parameters
clear
addpath('../src')  

%% Common parameters for both discretizations
%%% Frequency
omega = 2^4;
freq = omega/(2*pi); % Hz

%%% level of relative noise
delta_noise = 0.05;         

%%% domain parameters
a = 2;  % size of domain has to be an integer

% N_data = max(a*20+1,ceil(a*freq*8));
N_data = 163;
% N_data = 326;
h_data =  a/(N_data-1);

% N_solver = max(a*20+1,ceil(a*freq*8));
N_solver = 163;
% N_solver = 326;
h_solver =  a/(N_solver-1);

npml_data = max(20,round(2.5/(freq*h_data)));              % #pml layers
npml_solver = max(20,round(2.5/(freq*h_solver)));          % #pml layers

sigmaMax = 80;                            % intensity of the pml absorbtion

order_data = 4;                    % order of data generator
order_solver = 4;                    % order of inverse solver

%%% optimization parameters

alpha_i = 0;                  % initial guess
dTol = 1e-5;                  % first order optimality
maxIter = 1000;

%%% medium parameters
delta_m = 0.5;
rn = 0.3*sqrt(2);


%%% boundary parameters
xb = 0.0; yb = 0.0; rb = a/2;

%%% Print some parameters

lbd = 1/(sqrt(1+delta_m)*freq);

fprintf('freq = %.2f Hz\n', freq)
fprintf('lbd  = %.2e \n', lbd)
fprintf('h_data = %.2e \n', h_data)
fprintf('h_solver = %.2e \n', h_solver)

fprintf('N_data = %d\n', N_data)
fprintf('N_solver = %d\n', N_solver)

fprintf('PPW_data = %d\n', floor(lbd/h_data))
fprintf('PPW_solver = %d\n', floor(lbd/h_solver))

fprintf('npml_data = %d\n', npml_data)
fprintf('npml_solver = %d\n', npml_solver)

fprintf('order_data = %d\n', order_data)
fprintf('order_solver = %d\n', order_solver)

fprintf('Noise = %.2f\n',delta_noise)


%%% Source and receiver positions
% Number of directions
% Ntheta_i = max(80,floor(2*omega));
% Ntheta_r = max(80,floor(2*omega));
Ntheta_i = 180;
Ntheta_r = 180;


dtheta_i = 2*pi/(Ntheta_i);
dtheta_r = 2*pi/(Ntheta_r);

% incident directions
theta_i = linspace(pi, 3*pi-dtheta_i, Ntheta_i);
d = [cos(theta_i).' sin(theta_i).'];

% sampling the wavefield in the interpolated data points
theta_r = linspace(0, 2*pi-dtheta_r, Ntheta_r);
r = [cos(theta_r).' sin(theta_r).'];

points_query = [xb yb] + rb*r;


fprintf('Ni = %d \n', Ntheta_i)
fprintf('Nr = %d \n', Ntheta_r)

%% Generate data
%%% setup the model and the domain

% size of the model in interior domain
nxi  = N_data;
nyi  = N_data;
% ni   = nxi*nyi;

xi = h_data*(0:nxi-1) - 0.5*a;
yi = h_data*(0:nyi-1) - 0.5*a;

% [Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
nx = nxi + 2*npml_data;
ny = nyi + 2*npml_data;
% n  = nx*ny;

x  = [xi(1)+(-npml_data:-1)*h_data,xi,xi(end)+(1:npml_data)*h_data];
y  = [yi(1)+(-npml_data:-1)*h_data,yi,yi(end)+(1:npml_data)*h_data];

[X,Y] = meshgrid(x,y);

%%%  We generate the medium %%%

% we generate the Shepp-Logan phantom rescaled
xn = rn/sqrt(2);

Nn = length(xi(xi>-xn & xi<xn));

P = phantom('Modified Shepp-Logan',Nn);
P = delta_m/max(abs(P),[],'all')*P;

eta = zeros(nyi, nxi);

eta(yi>-xn & yi<xn, xi>-xn & xi<xn) = P;

eta = eta(:);

m = 1 + eta;

% extend the model to the simulation domain
eta_ext=ExtendModel(eta,nxi,nyi,npml_data);
mext=ExtendModel(m,nxi,nyi,npml_data);

figure(1); clf();
DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
% shading(gca,'interp');
title('Velocity');



% Helmholtz matrix without auxiliaries and with analytic 
%  differentiation of the pml profile
H1=HelmholtzMatrix(mext,nx,ny,npml_data,h_data,...
                   sigmaMax,order_data,omega,'compact_explicit');
H1 = H1';                    % The Helmholtz code produces adjoint matrix




% building the right hand sides
U_in =  exp(1i*omega*(X(:)*d(:,1).'+ Y(:)*d(:,2).'));
S    = bsxfun(@times, -omega^2*eta_ext, U_in);

% solving the equation
tic;
U = H1\S;
t_f = toc;

% printing the time of the solution
fprintf('Time elapsed of the computation = %.4e [s]\n',t_f );


%     % instead of interpolating at each iteration we compute the
%     % interpolation matrix and then this is reduced to a matrix vector
%     % multiplication
Projection_mat = zeros(Ntheta_r, nx, ny);
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

% this is our "real data"
scatter = reshape(Projection_mat, Ntheta_r, nx*ny)*U;


%% Save/load the data
% save(fullfile('..','data',...
%     ['IS_scatter_FD_SheppLogan_order',num2str(order_data),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_N',num2str(N_data),'_a',num2str(a),...
%         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.mat']),...
%         'eta','scatter','t_f');
%     
% load(fullfile('..','data',...
%     ['IS_scatter_FD_SheppLogan_order',num2str(order_data),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_N',num2str(N_data),'_a',num2str(a),...
%         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.mat']),...
%         'scatter');


%% add noise to the data
noise = delta_noise*scatter.*(2*binornd(1,0.5,size(scatter))-1);
scatter = scatter + noise;


%% Reconstruct medium
%%% setup the model and the domain

% size of the model in interior domain
nxi  = N_solver;
nyi  = N_solver;
% ni   = nxi*nyi;

xi = h_solver*(0:nxi-1) - 0.5*a;
yi = h_solver*(0:nyi-1) - 0.5*a;

% [Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
nx = nxi + 2*npml_solver;
ny = nyi + 2*npml_solver;
% n  = nx*ny;

x  = [xi(1)+(-npml_solver:-1)*h_solver,xi,xi(end)+(1:npml_solver)*h_solver];
y  = [yi(1)+(-npml_solver:-1)*h_solver,yi,yi(end)+(1:npml_solver)*h_solver];

[X,Y] = meshgrid(x,y);


%  we create a structure (this makes our lives easier)
% this should be encapsulated but for now we just use it.

% putting together all the different properties
properties.nx = nx;
properties.ny = ny;

properties.nxi = nxi;
properties.nyi = nyi;

properties.npml = npml_solver;
properties.omega = omega;

properties.Ntheta_i = Ntheta_i;
properties.Ntheta_r = Ntheta_r;

% intensity of the pml absorbtion
properties.sigmaMax = sigmaMax;
properties.h = h_solver;

% order of accuracy
% we use a different order to avoid the inverse crime
properties.order = order_solver;


%%%  We generate the medium %%%

% we generate the Shepp-Logan phantom rescaled
xn = rn/sqrt(2);

Nn = length(xi(xi>-xn & xi<xn));

P = phantom('Modified Shepp-Logan',Nn);
P = delta_m/max(abs(P),[],'all')*P;

eta = zeros(nyi, nxi);

eta(yi>-xn & yi<xn, xi>-xn & xi<xn) = P;

eta = eta(:);

m = 1 + eta;

% extend the model to the simulation domain
eta_ext=ExtendModel(eta,nxi,nyi,npml_solver);
mext=ExtendModel(m,nxi,nyi,npml_solver);

figure(1); clf();
DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
% shading(gca,'interp');
title('Velocity');

% building the incident plane wave
U_in =  exp(1i*omega*(X(:)*d(:,1).'+ Y(:)*d(:,2).'));
properties.U_in = U_in;

Projection_mat = zeros(Ntheta_r, nx, ny);
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
properties.project_mat = sparse(reshape(Projection_mat, Ntheta_r, nx*ny));








% define the misfit function (it provides both the misfit and the 
% gradient
J = @(x) misfit_FD(scatter, x, properties);

% prepare the options for the optimization loop 
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations', maxIter,...
    'OptimalityTolerance', dTol, ...
    'Display', 'iter-detailed');


% running the optimization routine
tic
[result,fval,exitflag,output] = fminunc(J,alpha_i*eta,options);
t_n = toc;

% printing the time of the solution
fprintf('Time elapsed of the optimization = %.2f [s]\n', t_n );


% plotting the result
figure(6);
clf();
subplot(1,3,1)
imagesc(xi,yi,reshape(eta, nxi, nyi)); pbaspect([1 1 1]); 
set(gca,'YDir','normal'); colorbar(); % shading interp;
title('Exact')
subplot(1,3,2)
imagesc(xi,yi,reshape(result, nxi, nyi)); pbaspect([1 1 1]);
set(gca,'YDir','normal'); colorbar(); % shading interp;
title('Reconstruction')
subplot(1,3,3)
imagesc(xi,yi,reshape(result-eta, nxi, nyi)); pbaspect([1 1 1]);
set(gca,'YDir','normal'); colorbar(); % shading interp;
title('Error')



% print and save
% set(gcf, 'Position',  [100, 100, 1500, 400])
% 
% % saveas(gcf,fullfile('..','plots',...
% %     ['IS_n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
% %     '_An',num2str(delta_m),'_rn',num2str(rn),...
% %         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
% %         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.eps']),'epsc')
% %     
% % print(fullfile('..','plots',...
% %     ['IS_n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
% %     '_An',num2str(delta_m),'_rn',num2str(rn),...
% %         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
% %         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.pdf']),'-dpdf');
% 
% saveas(gcf,fullfile('..','plots',...
%     ['IS_n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
%         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.fig']));
%               
% save(fullfile('..','data',...
%     ['IS_n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),....
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
%         '_Nthetai-r',sprintf('-%i',[Ntheta_i,Ntheta_r]),'.mat']),...
%         'eta','result','scatter','output','t_n');
