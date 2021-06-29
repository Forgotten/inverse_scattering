%% setup scaling parameters
clear
addpath('../src')  

%% butterfly and helmholtz parameters
nn = 2;

    freq = 10; % Hz
    omega = 2*pi*freq;

    delta_m = -0.2; % amplitude
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

    %%
    % normpdf with amplitude scale = 1 fixed
    normal = @(x, mu, sigma)(normpdf(x, mu, sigma)*sqrt(2*pi)*sigma);

    % number of scatterers (uniform between [2,4])

    eta = (0.2*normpdf(Xi(:),0.5,5*h).*normpdf(Yi(:),0.5,5*h));
    centres = 0.4*rand(nn,2) -0.2 ;
    
    sigma_h = 2*h;
    
    for ii = 1:nn
        if ii == 1
            eta = delta_m*normal(Xi(:),centres(ii,1),sigma_h).*normal(Yi(:),centres(ii,2),sigma_h);
        else
            eta = eta+delta_m*normal(Xi(:),centres(ii,1),sigma_h).*normal(Yi(:),centres(ii,2),sigma_h);
        end
    end
    % extend the model to the simulation domain
    m = 1 + eta;
    eta_ext=ExtendModel(eta,nxi,nyi,npml);
    mext=ExtendModel(m,nxi,nyi,npml);

    figure(1)
    DisplayField(1./sqrt(m),xi,yi);
    title('Velocity');

    % order of accuracy
    order = 8;

    % intensity of the pml absorbtion (the profile is quadratic from 0 to
    % sigmaMax
    sigmaMax = 80;

    %% Helmholtz matrix without auxiliaries and with analytic differentiation of the pml profile
    tic;
    H1=HelmholtzMatrix(mext,nx,ny,npml,h,sigmaMax,order,omega,'compact_explicit');
    t1=toc;
    [L1 U1 P1 Q1]=lu(H1);
    t2=toc;

    foo=whos('L1');
    fprintf('Time to compute the LU factorization for the compact analytic formulation = %.2fs\n',t2-t1);
    fprintf('Size of the L factor = %.1fMb\n',foo.bytes/1024^2);

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
    U = Q1*(U1\(L1\(P1*S)));

    % sampling the wavefield in the interpolated data points

    scatter = zeros(Ntheta, Ntheta);
    theta_r = linspace(0, 2*pi-dtheta, Ntheta);
%     theta_r = theta
    r = [cos(theta_r).' sin(theta_r).'];
    
    points_query = 0.5*r;
    
    for ii= 1:Ntheta
        scatter(:, ii) = interp2(x,y,reshape(U(:,ii), nx, ny), points_query(:,1), points_query(:,2));
    end

    %% reconstruction
    [pa,qa] = meshgrid(theta, theta);
    p1 = cos(pa);    p2 = sin(pa);
    q1 = cos(qa);    q2 = sin(qa);
    p1 = p1(:);    p2 = p2(:);
    q1 = q1(:);    q2 = q2(:);
    
    F = (-omega^2*exp(1i*omega*0.5)/sqrt(0.5))* (exp(1i*omega*( (p1 - q1).*X(:)'+...
                        (p2 - q2).*Y(:)')) );
    F_eta= F*eta_ext;
    
    figure(2); clf;
    subplot(1,2,1);
    imagesc(reshape(real(F_eta), 80, 80))
    subplot(1,2,2);
    imagesc(reshape(real(scatter), 80, 80))
    
    eps_reg = 1.e-3;
    filter = @(x) F'*(F*x) + eps_reg*x;
    
    a = gmres(filter, F'*scatter(:), 10,1e-6,10);
    
    figure(3); clf();
    subplot(1,3,1);
    imagesc(reshape(real(F'*scatter(:)), 120, 120))
    subplot(1,3,2);
    imagesc(reshape(real(a), 120, 120))
    subplot(1,3,3);
    imagesc(reshape(real(eta_ext), 120, 120))
 

