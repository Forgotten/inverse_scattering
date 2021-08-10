function [mis,  dmis] =  misfit_multifreq(Data_array, eta, props)

    % unrolling the data contained in the structure 
    nxi = props.nxi;
    nyi = props.nyi;
    
    nx = props.nx;
    ny = props.ny;
    
    npml = props.npml;
    % in this case we need an array of omega
    omega_array = props.omega_array;
    
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
    eta_ext=ExtendModel(eta,nxi,nyi,npml);
    m_ext=ExtendModel(m,nxi,nyi,npml);
    
    
    % variable to accumulate the misfit
    mis = 0;
    
    if nargout > 1
        dmis = 0;
    end
    
    for ii = 1:length(omega_array)
        
        omega = omega_array(ii);
        Data  = squeeze(Data_array(ii,:,:));
        
        %% Helmholtz matrix without auxiliaries and with analytic differentiation of the pml profile
        H1=HelmholtzMatrix(m_ext,nx,ny,npml,h,sigmaMax,...
                           order,omega,'compact_explicit');

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
        U = H1\S; 

        % sampling the wavefield in the interpolated data points
        scatter = Projection_mat*U;

        residual = Data - scatter;

        mis_omega =  0.5*sum(abs(residual).^2,'all');
        mis = mis + mis_omega;

        if nargout > 1

            % computing the total wavefield  
            U_total = U + U_in;

            % computing the rhs for the adjoint system
            rhs_adj = omega^2*(Projection_mat.')*residual;

            % solving the adjoint system
            W_adj = (H1')\rhs_adj;

            % computing the gradient
            grad = real(sum(conj(U_total).*W_adj, 2));

            % reshaping and extrating the gradient
            grad = reshape(grad, nx,ny);
            grad = grad(npml+1:npml+nxi, npml+1:npml+nyi);

            % computing the gradient for a given frequency
            d_mis_omega = grad(:);

            dmis = dmis + d_mis_omega;
        end 

    end 
    
end
    