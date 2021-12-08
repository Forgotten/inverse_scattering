% This code use the Lippman-Schwinger solver for data generation and 
% finite difference for stochastic gradient descent

%% setup path
clear
addpath('../src')  

%% Common parameters for both methods
%%% Frequency
freq = 2^4; % Hz            % wavelength = 1/freq
omega = 2*pi*freq;

delta_noise = 0.0;         % level of relative noise

%%% domain parameters
a = 2;                        % domain size

N_LS = a*freq*2^4;            % #grid points of Lippman-Schwinger
h_LS = a/N_LS;                % step size

N_FD = a*freq*32;              % #grid points of finite difference
h_FD = a/(N_FD-1);             


%%% optimization parameters
order = 4;                    % order of FD scheme

%%% medium parameters
nn = 2;                              % number of bumps
MDName = [0,0];                      % shapes of supports (polygon or disk or star polygon)
% radius of support
rn = [0.1,0.2];
% center of heterogeneity
centres = [-0.2,-0.2;0.1,0.1];       % 1st/2nd column: x/y coordinate
% center of heterogeneity
delta_m_nn = [0.5,0.3];              % amplitude
delta_m = max(delta_m_nn); 


%%% Print some parameters

lbd = 1/(sqrt(1+delta_m)*freq);      % lower bound of wavelength

fprintf('freq = %.2f Hz\n', freq)
fprintf('lbd  = %.2e \n', lbd)
fprintf('h_LS = %.2e \n', h_LS)
fprintf('h_FD = %.2e \n', h_FD)

fprintf('N_LS = %d\n', N_LS)
fprintf('N_FD = %d\n', N_FD)
fprintf('PPW_LS = %d\n', floor(lbd/h_LS))
fprintf('PPW_FD = %d\n', floor(lbd/h_FD))


%%% boundary parameters
xb = 0.0; yb = 0.0; rb = 0.4;


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
Ntheta_s = 80;
Ntheta_r = 80;
% Number of directions
Ntheta_i = 40;
Ntheta_o = 40;

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







%% Generate data
%%% Lippmann-Schwinger parameters

% domain parameters
x_LS = -a/2:h_LS:a/2-h_LS;
y_LS = -a/2:h_LS:a/2-h_LS;

X_LS = repmat(x_LS',1,N_LS);
Y_LS = repmat(y_LS,N_LS,1);

% medium matrix
eta_LS = zeros(N_LS,N_LS);
for ii = 1:nn
    eta_LS = eta_LS + delta_m_nn(ii)*...
        polygon_bump(X_LS,Y_LS,centres(ii,1),centres(ii,2),rn(ii),MDName(ii));
end
m = 1 + eta_LS;

figure(1); clf();
DisplayField(1./sqrt(m.'),x_LS,y_LS); set(gca,'YDir','normal');
title('Velocity');

tic;
%%% Lippmann-Schwinger operator
LS = LippmannSchwinger_precompute(x_LS,y_LS,omega,eta_LS);


%%% solving the equation
U = zeros(N_LS*N_LS,Ntheta_s*Ntheta_i);
for ii = 1:Ntheta_s*Ntheta_i
    
    %%% Generate source
    S_LS = Sfunc(X_LS(:)-pos(ii,1) , Y_LS(:)-pos(ii,2) , dir(ii,1) , dir(ii,2) );
    
    
    %%% Plot source
%     S_plot = reshape(S_LS,N_LS,N_LS);
%     figure(2); clf();
%     DisplayField(S_plot.',x_LS,y_LS); set(gca,'YDir','normal');
%     title('Source'); pause;

    %%% Building the incident wave
    u_inc = LS.apply_Green(S_LS);
    
    
    %%% building the right hand-side
    rhsDual = -omega^2*eta_LS.*u_inc;
    
    
    %%% solving the Lippmann-Schwinger equation
    sigma = LS\rhsDual(:);
    
    
    %%% computing the wavefield
    u_sca = LS.apply_Green(sigma);
    u_tot = u_inc + u_sca;
    
    U(:,ii) = u_tot(:);
    
    
    %%% Plot solution
%     figure(3); clf();
%     DisplayField(u_tot.',x_LS,y_LS); set(gca,'YDir','normal');
%     title('Solution'); pause;

end
t_f = toc;

% printing the time of the solution
fprintf('Time elapsed of the computation = %.4e [s]\n',t_f );


%%% Perform Husimi transform
husimi_mat = phi( pos_r(:,1)-X_LS(:).' , pos_r(:,2)-Y_LS(:).' , dir_r(:,1) , dir_r(:,2) );


%%% this is our "real data" %%%
scatter = abs(h_LS^2*husimi_mat*U).^2;


%% Save/load the data
save(fullfile('..','data',...
    ['scatter_LS_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
        '_f',num2str(freq),'_N',num2str(N_LS),'_a',num2str(a),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
        'eta_LS','scatter');
    
% load(fullfile('..','data',...
%     ['scatter_LS_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
%         '_f',num2str(freq),'_N',num2str(N_LS),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
%         'scatter');

%%% add noise to the data
noise = delta_noise*scatter.*(2*binornd(1,0.5,size(scatter))-1);
scatter = scatter + noise;




%% Reconstruct medium
%%% setup the model and the domain

% size of the model in interior domain
nxi  = N_FD;
nyi  = N_FD;
ni   = nxi*nyi;

xi = h_FD*(0:nxi-1) - 0.5*a;
yi = h_FD*(0:nyi-1) - 0.5*a;

[Xi,Yi] = meshgrid(xi,yi);

% size of the simulation domain
npml = 20;
nx = nxi + 2*npml;
ny = nyi + 2*npml;
% n  = nx*ny;

x  = [xi(1)+(-npml:-1)*h_FD,xi,xi(end)+(1:npml)*h_FD];
y  = [yi(1)+(-npml:-1)*h_FD,yi,yi(end)+(1:npml)*h_FD];

[X,Y] = meshgrid(x,y);

% intensity of the pml absorbtion
sigmaMax = 80;



%%%  we create a structure (this makes our lives easier)
% this should be encapsulated but for now we just use it.

% putting together all the different properties
properties.nx = nx;
properties.ny = ny;
% properties.n = n;

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
properties.h = h_FD;

% order of accuracy
properties.order = order;



%%%  We generate the true medium 
eta = zeros(ni,1);
for ii = 1:nn
    eta = eta + delta_m_nn(ii)*...
        polygon_bump(Xi(:),Yi(:),centres(ii,1),centres(ii,2),rn(ii),MDName(ii));
end
% extend the model to the simulation domain
m = 1 + eta;

figure(4); clf();
DisplayField(1./sqrt(m),xi,yi); set(gca,'YDir','normal');
title('Velocity');
  


%%% computing source for probing

% building the source matrix
S = Sfunc(X(:)-pos(:,1).' , Y(:)-pos(:,2).' , dir(:,1).' , dir(:,2).' );

S_plot = reshape(S(:,1),ny,nx); 
S_plot = S_plot(npml+1:npml+nxi, npml+1:npml+nyi);
figure(5); clf();
DisplayField(S_plot,xi,yi); set(gca,'YDir','normal');
title('Source');

% we save the source matrix into the properties structure
properties.S = S;

%%% Compute the Husimi matrix
husimi_mat = phi( pos_r(:,1)-X(:).' , pos_r(:,2)-Y(:).' , dir_r(:,1) , dir_r(:,2) );

% we save the husimi matrix into the properties structure
properties.husimi_mat = husimi_mat;



%%% SGD
% define the misfit function (it provides the gradient)
% J = @(x,reduction,rnd_num) misfit_husimi_sgd(scatter, x, properties, reduction, rnd_num);

iter_max = 200;

batch_size = round(0.01*Ntheta_s*Ntheta_i);

output_step = 50;

% learning rate scheduler
lr = 100;
lr_step = 200;
alpha_lr = 0.5;

% initial guess
alpha_i = 0;                  % initial guess
eta_inv = alpha_i*eta;

% print initial loss and learning rate
loss = misfit_husimi_sgd(scatter, eta_inv, properties, 'mean');
fprintf('Step: %d \t  Loss: %.8f \t  Learning Rate: %.5f \n',...
    0, loss, lr)

% SGD
rng(0);
ii = 1;
while ii <= iter_max
    % random source position/direction
    gamma = randperm(Ntheta_s*Ntheta_i,batch_size);
    % compute loss and gradient
    [loss, dloss] = misfit_husimi_sgd(scatter, eta_inv, properties, 'sum', gamma);
    % update
    eta_inv = eta_inv - lr/batch_size*dloss;
    
    % print
    if rem(ii,output_step) == 0
        loss = misfit_husimi_sgd(scatter, eta_inv, properties, 'mean');
        fprintf('Step: %d \t  Loss: %.8f \t  Learning Rate: %.5f \n',...
            ii, loss, lr)
    end
    
    ii = ii + 1;
    % update learning rate
    if rem(ii,lr_step) == 1
        lr = lr*alpha_lr;
    end
end
fprintf('Done! \n')

result = eta_inv;



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



%% print and save
set(gcf, 'Position',  [100, 100, 1500, 400])

saveas(gcf,fullfile('..','plots',...
    ['n_SGD_order',int2str(order),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),...
    '_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
        '_f',num2str(freq),'_Nfd-ls',sprintf('-%i',[N_FD,N_LS]),'_a',num2str(a),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.eps']),'epsc')
    
% print(fullfile('..','plots',...
%     ['n_SGD_order',int2str(order),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),...
%     '_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
%         '_f',num2str(freq),'_Nfd-ls',sprintf('-%i',[N_FD,N_LS]),'_a',num2str(a),...
%         '_sigma',int2str(-log2(sigma0)),...
%         '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.pdf']),'-dpdf');

savefig(fullfile('..','plots',...
    ['n_SGD_order',int2str(order),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),...
    '_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
        '_f',num2str(freq),'_Nfd-ls',sprintf('-%i',[N_FD,N_LS]),'_a',num2str(a),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.fig']));
              
save(fullfile('..','data',...
    ['n_SGD_order',int2str(order),'_etai',num2str(alpha_i),'_dNoise',num2str(delta_noise),...
    '_medium',sprintf('-%i',MDName),'_An',sprintf('-%.2f',delta_m_nn),'_rn',sprintf('-%.2f',rn),...
        '_f',num2str(freq),'_Nfd-ls',sprintf('-%i',[N_FD,N_LS]),'_a',num2str(a),...
        '_sigma',int2str(-log2(sigma0)),...
        '_Nthetai-s-o-r',sprintf('-%i',[Ntheta_i,Ntheta_s,Ntheta_o,Ntheta_r]),'.mat']),...
        'eta','result','scatter');

