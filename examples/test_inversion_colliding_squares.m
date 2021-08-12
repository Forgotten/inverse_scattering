%% setup scaling parameters
clear
addpath('../src')  

% Prototype for full-waveform inversion with multiple 
% frequencies (work in progress)

sigma = 0.01;
L = 4;
s = 5; % make leaf node several wavelengths?
N = (2^L)*s;
h =  1/(N-1);

%% setup the model and the domain

% background wavespeed
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
properties.nx = nx;      properties.ny = ny;
properties.nxi = nxi;    properties.nyi = nyi;
properties.X = X;        properties.Y = Y;
properties.npml = npml;  properties.N = N;

% intensity of the pml absorbtion
properties.sigmaMax = sigmaMax;  
properties.h = h;

% order of accuracy
% we use a different order to avoid the inverse crime
properties.order = 4;

%%  We load the data from the experiments

folder = '../../Inverse_Scattering_ML_TF2/data/colliding_data_L4s5_multifreq_square_10_3_h_freq_2.5_5_10/';

eta_array = readmatrix(strcat(folder, 'eta_freq_10.0.csv'));

scatter_array = zeros(3, size(eta_array,1),  size(eta_array,2));             
               
scatter_array(1,:,:) = readmatrix(strcat(folder,'scatter_real_freq_2.5.csv')) + ...
                     1i*readmatrix(strcat(folder,'scatter_imag_freq_2.5.csv'));
scatter_array(2,:,:) = readmatrix(strcat(folder,'scatter_real_freq_5.0.csv')) + ...
                     1i*readmatrix(strcat(folder,'scatter_imag_freq_5.0.csv'));
scatter_array(3,:,:) = readmatrix(strcat(folder,'scatter_real_freq_10.0.csv')) + ...
                     1i*readmatrix(strcat(folder,'scatter_imag_freq_10.0.csv'));

% we choose one particular experiment of the series
idx = 9;

eta = reshape(eta_array(idx,:), nxi, nyi).';
eta = eta(:);


m = 1 + eta;
               
figure(1); clf();
subplot(1,2,1);
DisplayField(1./sqrt(m),xi,yi);
title('Velocity');
subplot(1,2,2);
DisplayField(eta,xi,yi);
title('Perturbation');

%% we build the projection matrix 

Ntheta = N;
dtheta = 2*pi/(Ntheta);

theta = linspace(pi, 3*pi-dtheta, Ntheta);
d = [cos(theta).' sin(theta).'];

theta_r = linspace(0, 2*pi-dtheta, Ntheta);
r = [cos(theta_r).' sin(theta_r).'];

points_query = 0.5*r;

project_mat = zeros(Ntheta, nx, ny);

for ii = 1:nx
    for jj = 1:ny
        mat_dummy = zeros(nx,ny);
        mat_dummy(ii,jj) = 1;
        project_mat(:,ii,jj) = interp2(x,y,...
                                       reshape(mat_dummy, nx, ny),...
                                       points_query(:,1),...
                                       points_query(:,2));
    end
end

% properly reshaping and making it sparse
project_mat = sparse(reshape(project_mat, Ntheta, nx*ny));

% we save the projection matrix into the properties structure
properties.project_mat = project_mat;

% selecting inversion frequencies
frequencies = [2.5 5 10];

% we start from a zero guess
eta_0 = 0*eta;

for i = 1:3 

    freq = frequencies(i);
    omega = 2*pi*freq;
    properties.omega = omega;
    
    
    scatter = reshape(scatter_array(i,idx,:), nxi, nyi).';
    
    % plotting the far field
    figure(2); clf();
    subplot(1,2,1);
    imagesc(real(reshape(scatter,nxi,nyi)));
    title('real part of the far field');
    subplot(1,2,2);
    imagesc(imag(reshape(scatter,nxi,nyi)));
    title('imaginary part of the far field');
    
    %% running the optimization
    
    % define the misfit function (it provides both the misfit and the
    % gradient
    J = @(x) misfit(scatter, x, properties);
    
    % prepare the options for the optimization loop
    options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'SpecifyObjectiveGradient',true,...
        'MaxIterations', 100,...
        'OptimalityTolerance', 1e-5, ...
        'Display', 'iter-detailed');
    
    tic;
    % running the optimization routine
    [result,fval,exitflag,output] = fminunc(J,eta_0,options);
    t_f = toc;
    
    % printing the time of the solution
    fprintf('Time elapsed for the optimization = %.4e [s]\n',t_f );
    
    % plotting the result
    figure(5+i);
    clf();
    subplot(1,3,1)
    imagesc(reshape(eta, nxi, nyi))
    colorbar();
    title('Exact')
    subplot(1,3,2)
    imagesc(reshape(result, nxi, nyi))
    colorbar();
    title('Reconstruction')
    subplot(1,3,3)
    imagesc(reshape(result-eta, nxi, nyi))
    colorbar();
    title('Error')
    
    % we update the guess for each frequency
    eta_0 = result;
    
end
