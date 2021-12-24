
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
xb = 0.0; yb = 0.0; rb = 0.4;


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
fprintf('npml_data = %d\n', npml_solver)

fprintf('order_data = %d\n', order_data)
fprintf('order_solver = %d\n', order_solver)



%%% Source parameters
% Width
sigma0 = 2^(-2);
Cd = sqrt(2)*(sigma0/sqrt(pi))^(3/2);
% Prototype function
Sfunc = @(x,y,thetax,thetay) omega^(5/2)*Cd*...
            exp(-omega^2*sigma0^2/2*(x.^2+y.^2) + 1i*omega*(thetax.*x+thetay.*y));

        
%%% Gabor filter for Husimi transform
phi = @(x,y,thetax,thetay) 0.5*(omega/pi)^(3/2)*...
                exp(-omega/2*(x.^2+y.^2) + 1i*omega*(thetax.*x+thetay.*y));

            
%%% Source and receiver positions
% Number of sources/receivers
Ntheta_s = 96; 
Ntheta_r = 96; 
% Number of directions
Ntheta_i = 48; 
Ntheta_o = 48; 

% Ntheta_s = max(80,floor(1.5*omega));
% Ntheta_r = max(80,floor(1.5*omega));
% Ntheta_i = max(40,floor(0.75*omega));
% Ntheta_o = max(40,floor(0.75*omega));

% source position/direction
dtheta_s = 2*pi/Ntheta_s;
dtheta_i = pi/(Ntheta_i+1);

theta_s = repmat(linspace(pi, 3*pi-dtheta_s, Ntheta_s), Ntheta_i, 1);
theta_i = repmat(linspace(-pi/2+dtheta_i, pi/2-dtheta_i, Ntheta_i)', 1, Ntheta_s);
theta_i = theta_i + theta_s + pi;   

dir = [cos(theta_i(:)), sin(theta_i(:))];
pos = [xb+rb*cos(theta_s(:)) , yb+rb*sin(theta_s(:))];

% receiver position/direction
dtheta_r = 2*pi/Ntheta_r;
dtheta_o = pi/(Ntheta_o+1);

theta_r = repmat(linspace(0, 2*pi-dtheta_r, Ntheta_r), Ntheta_o, 1);
theta_o = repmat(linspace(-pi/2+dtheta_o, pi/2-dtheta_o, Ntheta_o)', 1, Ntheta_r);
theta_o = theta_o + theta_r;

dir_r = [cos(theta_o(:)), sin(theta_o(:))];
pos_r = [xb+rb*cos(theta_r(:)), yb+rb*sin(theta_r(:))];


fprintf('Ns = %d \n', Ntheta_s)
fprintf('Ni = %d \n', Ntheta_i)
fprintf('Nr = %d \n', Ntheta_r)
fprintf('No = %d \n', Ntheta_o)



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



% extend the model to the simulation domain
m = 1 + eta;
mext=ExtendModel(m,nxi,nyi,npml_data);

figure(1); clf();
DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
title('Velocity');

        
%%% computing source for probing

% building the source matrix
S = Sfunc(X(:)-pos(:,1).' , Y(:)-pos(:,2).' , dir(:,1).' , dir(:,2).' );


S_plot = reshape(S(:,1),ny,nx); 
S_plot = S_plot(npml_data+1:npml_data+nxi, npml_data+1:npml_data+nyi);
figure(2); clf();
DisplayField(S_plot,xi,yi); set(gca,'YDir','normal');
title('Source');


%%% Generate the Helmholtz matrix %%%

% Helmholtz matrix without auxiliaries and with analytic 
%  differentiation of the pml profile
H1=HelmholtzMatrix(mext,nx,ny,npml_data,h_data,...
                   sigmaMax,order_data,omega,'compact_explicit');
H1 = H1';                    % The Helmholtz code produces adjoint matrix

       

% solving the equation
tic;
U = H1\S;
t_f = toc;

% printing the time of the solution
fprintf('Time elapsed of the computation = %.2f [s]\n',t_f );

U_plot = reshape(U(:,1),ny,nx); 
U_plot = U_plot(npml_data+1:npml_data+nxi, npml_data+1:npml_data+nyi);
figure(3); clf();
DisplayField(U_plot,xi,yi); set(gca,'YDir','normal');
title('Solution');

%%% Perform Husimi transform %%%
husimi_mat = phi( pos_r(:,1)-X(:).' , pos_r(:,2)-Y(:).' , dir_r(:,1) , dir_r(:,2) );

%%% this is our "real data" 
scatter = abs(h_data^2*husimi_mat*U).^2;

%% Save/load the data
% save(fullfile('..','data',...
%     ['scatter_FD_SheppLogan_order',num2str(order_data),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_N',num2str(N_data),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
%         'eta','scatter','t_f');
%     
% load(fullfile('..','data',...
%     ['scatter_FD_SheppLogan_order',num2str(order_data),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_N',num2str(N_data),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
%         'scatter');

%% add noise to the data
noise = delta_noise*scatter.*(2*binornd(1,0.5,size(scatter))-1);
scatter = scatter + noise;


%% Reconstruct medium
%%% setup the model and the domain

% size of the model in interior domain
nxi  = N_solver;
nyi  = N_solver;
ni   = nxi*nyi;

xi = h_solver*(0:nxi-1) - 0.5*a;
yi = h_solver*(0:nyi-1) - 0.5*a;

[Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
nx = nxi + 2*npml_solver;
ny = nyi + 2*npml_solver;
% n  = nx*ny;

x  = [xi(1)+(-npml_solver:-1)*h_solver,xi,xi(end)+(1:npml_solver)*h_solver];
y  = [yi(1)+(-npml_solver:-1)*h_solver,yi,yi(end)+(1:npml_solver)*h_solver];

[X,Y] = meshgrid(x,y);


%%%  we create a structure (this makes our lives easier)
% this should be encapsulated but for now we just use it.

% putting together all the different properties
properties.nx = nx;
properties.ny = ny;
% properties.n = n;

properties.nxi = nxi;
properties.nyi = nyi;

properties.npml = npml_solver;
properties.omega = omega;

properties.Ntheta_s = Ntheta_s;
properties.Ntheta_i = Ntheta_i;
properties.Ntheta_r = Ntheta_r;
properties.Ntheta_o = Ntheta_o;

% intensity of the pml absorbtion
properties.sigmaMax = sigmaMax;
properties.h = h_solver;

% order of accuracy
properties.order = order_solver;

%%%  We generate the medium %%%

xn = rn/sqrt(2);

Nn = length(xi(xi>-xn & xi<xn));

P = phantom('Modified Shepp-Logan',Nn);
P = delta_m/max(abs(P),[],'all')*P;

eta = zeros(nyi, nxi);

eta(xi>-xn & xi<xn, xi>-xn & xi<xn) = P;

eta = eta(:);

% plot the model
% m = 1 + eta;
% 
% figure(4); clf();
% DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
% title('Velocity');


%%% computing source for probing

% building the source matrix
S = Sfunc(X(:)-pos(:,1).' , Y(:)-pos(:,2).' , dir(:,1).' , dir(:,2).' );

% S_plot = reshape(S(:,1),ny,nx); 
% S_plot = S_plot(npml_solver+1:npml_solver+nxi, npml_solver+1:npml_solver+nyi);
% figure(5); clf();
% DisplayField(S_plot,xi,yi); set(gca,'YDir','normal');
% title('Source');

% we save the source matrix into the properties structure
properties.S = S;

%%% Perform Husimi transform
husimi_mat = phi( pos_r(:,1)-X(:).' , pos_r(:,2)-Y(:).' , dir_r(:,1) , dir_r(:,2) );

% we save the husimi matrix into the properties structure
properties.husimi_mat = husimi_mat;




%%% Nonlinear least square
% define the misfit function (it provides both the misfit and the 
% gradient

% J = @(x) misfit_husimi_prototype(scatter, x, properties);
J = @(x) misfit_husimi(scatter, x, properties);


% prepare the options for the optimization loop 
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations', maxIter,...
    'OptimalityTolerance', dTol, ...
    'Display', 'iter-detailed');


tic
% running the optimization routine
[result,fval,exitflag,output] = fminunc(J,alpha_i*eta,options);
t_n = toc;

% printing the time of the solution
fprintf('Time elapsed of the optimization = %.2f [s]\n', t_n );


% plotting the result
figure(6);
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



% print and save
% set(gcf, 'Position',  [100, 100, 1500, 400])
% 
% % saveas(gcf,fullfile('..','plots',...
% %     ['n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
% %     '_An',num2str(delta_m),'_rn',num2str(rn),...
% %         '_f',num2str(freq),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
% %         '_sigma',int2str(-log2(sigma0)),...
% %         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.eps']),'epsc')
% %     
% % print(fullfile('..','plots',...
% %     ['n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
% %     '_An',num2str(delta_m),'_rn',num2str(rn),...
% %         '_f',num2str(freq),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
% %         '_sigma',int2str(-log2(sigma0)),...
% %         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.pdf']),'-dpdf');
% 
% saveas(gcf,fullfile('..','plots',...
%     ['n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),...
%     '_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.fig']));
%               
% save(fullfile('..','data',...
%     ['n_SheppLogan_order',sprintf('-%i',[order_data,order_solver]),...
%     '_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),'_dTol',int2str(-log10(dTol)),...
%     '_An',num2str(delta_m),'_rn',num2str(rn),...
%         '_k',num2str(omega),'_Ndata-solver',sprintf('-%i',[N_data,N_solver]),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
%         'eta','result','scatter','output','t_n');

