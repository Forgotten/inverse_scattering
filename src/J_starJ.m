function jj = J_starJ(eta, delta_eta, props, epsilon)
% eta is the backgroun media
% delta_eta is the perturbation

if nargin < 4
    epsilon = 0;
end

% unrolling the data contained in the structure
nxi = props.nxi;
nyi = props.nyi;

nx = props.nx;
ny = props.ny;

npml = props.npml;
omega = props.omega;

N = props.N;

X = props.X;
Y = props.Y;

Projection_mat = props.project_mat;

% intensity of the pml absorbtion (the profile is quadratic from 0 to
% sigmaMax
sigmaMax = props.sigmaMax;
h = props.h;

% order of accuracy
order = props.order;


%% preparing the media
% extend the model to the simulation domain
m = 1 + eta;

eta_ext = ExtendModel(eta,nxi,nyi,npml);
m_ext = ExtendModel(m,nxi,nyi,npml);

% we extend the perturbation too
delta_eta_ext = ExtendModel(delta_eta,nxi,nyi,npml);

%% Helmholtz matrix without auxiliaries and with analytic 
% differentiation of the pml profile
H = HelmholtzMatrix(m_ext,nx,ny,npml,h,sigmaMax,...
                    order,omega,'compact_explicit');

%% computing the 

% computing the number of directions
Ntheta = N;
dtheta = 2*pi/(Ntheta);

theta = linspace(pi, 3*pi-dtheta, Ntheta);
d = [cos(theta).' sin(theta).'];

% building the right hand sides
U_in =  exp(1i*omega*(X(:)*d(:,1).'+ Y(:)*d(:,2).'));
S    = bsxfun(@times, -omega^2*eta_ext, U_in);

% solving the equation
U = H\S;

% computing the total wavefield
U_total = U + U_in;

% building the right hand side for the Born approximation
S_born = bsxfun(@times, -omega^2*delta_eta_ext, U_total);

% solving the equation
U_born = H\S_born;

% sampling the wavefield in the interpolated data points
scatter = Projection_mat*U_born;

% computing the rhs for the adjoint system
rhs_adj = -omega^2*(Projection_mat.')*scatter;

% solving the adjoint system
W_adj = (H')\rhs_adj;

% computing the gradient
grad = real(sum(conj(U_total).*W_adj, 2));

% reshaping and extrating the gradient
grad = reshape(grad, nx,ny);
grad = grad(npml+1:npml+nxi, npml+1:npml+nyi);

% returnig the gradient
jj = grad(:) + epsilon*delta_eta;




end
