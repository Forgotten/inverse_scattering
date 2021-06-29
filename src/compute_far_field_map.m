%% setup the model and the domain
clear
% background wavespeed
N = 100;
c = ones(N,N);

h =  1/(N-1) ;
% size of the model
nxi  = size(c,2);
nyi  = size(c,1);
ni   = nxi*nyi;

npml = 20;
% size of the simulation domain
nx = nxi + 2*npml;
ny = nyi + 2*npml;
n  = nx*ny;

xi = h*(0:nxi-1) - 0.5; % domain [-0.5, 0.5]
yi = h*(0:nyi-1) - 0.5; % domain [-0.5, 0.5]
x  = [xi(1)+(-npml:-1)*h xi xi(end)+(1:npml)*h];
y  = [yi(1)+(-npml:-1)*h yi yi(end)+(1:npml)*h];

[X,Y]   = meshgrid(x,y);
[Xi,Yi] = meshgrid(xi,yi);

%%

% extend the model to the simulation domain
eta = (0.2*normpdf(Xi(:),0.1,5*h).*normpdf(Yi(:),0.1,5*h));
%eta = 0.2*exp((Xi(:)-0.5).^2/(3*h)^2).*exp((Yi(:)-0.5).^2/(3*h)^2);

m = 1 + eta;
eta_ext = ExtendModel(eta,nxi,nyi,npml);
mext = ExtendModel(m,nxi,nyi,npml);

figure(1)
DisplayField(1./sqrt(m),xi,yi);
title('Velocity');

% frequency = 5Hz
omega = 2*pi*nyi/10;

% order of accuracy
order = 8;

% intensity of the pml absorbtion (the profile is quadratic from 0 to
% sigmaMax
sigmaMax = 80;

%% source term (in the simulation domain)

% incoming direction
d_s = [0.0, 1.0];

u_in = exp(1i*omega*(d_s(1)*X(:)+ d_s(2)*Y(:)));
s = -omega^2*eta_ext.*u_in;

% plotting the source
figure(2)
DisplayField(s,x,y, npml);


%% Helmholtz matrix without auxiliaries and with analytic differentiation of the pml profile
tic;

H1=HelmholtzMatrix(mext,nx,ny,npml,h,sigmaMax,order,omega,'compact_explicit');
t1=toc;
[L1 U1 P1 Q1]=lu(H1);
t2=toc;

foo=whos('L1');
fprintf('Time to compute the LU factorization for the compact analytic formulation = %.2fs\n',t2-t1);
fprintf('Size of the L factor = %.1fMb\n',foo.bytes/1024^2);

%% solve
u1=Q1*(U1\(L1\(P1*s)));

%% display
figure(9)
DisplayField(u1,x,y,npml);
title('Compact/exact')

%% computing the farfield data

% computing the number of directions
% Ntheta = 1000;
Ntheta = 100;
dtheta = 2*pi/(Ntheta);

theta = linspace(pi, 3*pi-dtheta, Ntheta);
d_s = [cos(theta).' sin(theta).'];

% building the right hand sides
U_in =  exp(1i*omega*(X(:)*d_s(:,1).'+ Y(:)*d_s(:,2).'));
S    = bsxfun(@times, -omega^2*eta_ext, U_in);

% solving the equation
U =Q1*(U1\(L1\(P1*S)));

% sampling the wavefield in the interpolated data points

Lambda = zeros(Ntheta, Ntheta);
% Lambda_{r,s}

theta = dtheta*(0:Ntheta-1);
d_rcv = [cos(theta).' sin(theta).'];

% r
points_query = 0.5*d_rcv;

for ii= 1:Ntheta
    Lambda(:, ii) = interp2(x,y,reshape(U(:,ii), nx, ny), points_query(:,1), points_query(:,2));
end

figure(5); clf();
imagesc(real(Lambda))
title('real(scattering)')

figure(6); clf();
imagesc(imag(Lambda))
title('imag(scattering)')
