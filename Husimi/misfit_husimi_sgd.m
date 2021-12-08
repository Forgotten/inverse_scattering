function [mis,  dmis] =  misfit_husimi_sgd(Data, eta, props, reduction, rnd_num)
    
    % unrolling the data contained in the structure 
    nxi = props.nxi;
    nyi = props.nyi;
    
    nx = props.nx;
    ny = props.ny;
%     n = props.n;
    
    npml = props.npml;
    omega = props.omega;
    
    Ntheta_s = props.Ntheta_s;
    Ntheta_i = props.Ntheta_i;
    Ntheta_r = props.Ntheta_r;
    Ntheta_o = props.Ntheta_o;
    
    husimi_mat = props.husimi_mat;
    
    S = props.S;
    if nargin >= 5
        S = S(:,rnd_num);
        
        Data = Data(:,rnd_num);
    end
    
    % intensity of the pml absorbtion (the profile is quadratic from 0 to
    % sigmaMax
    sigmaMax = props.sigmaMax;
    h = props.h; 
    
    % order of accuracy
    order = props.order;
    
    %% preparing the media 
    % extend the model to the simulation domain
    m = 1 + eta;
%     eta_ext=ExtendModel(eta,nxi,nyi,npml);
    m_ext = ExtendModel(m,nxi,nyi,npml);

    %% Helmholtz matrix without auxiliaries and with analytic differentiation of the pml profile
    H1=HelmholtzMatrix(m_ext,nx,ny,npml,h,sigmaMax,...
                        order,omega,'compact_explicit');
    H1 = H1';

    %% computing the farfield data

    % solving the equation
    U = H1\S; 

    % sampling the wavefield in the interpolated data points
    phi_u = h^2*husimi_mat*U;
    scatter = abs(phi_u).^2;
    
    residual = Data - scatter;
    
    if reduction == "sum"
        mis =  0.5*sum((residual).^2,'all')/(Ntheta_r*Ntheta_o);
    elseif reduction == "mean"
        mis = 0.5*sum((residual).^2,'all')/(Ntheta_s*Ntheta_i*Ntheta_r*Ntheta_o);
    else
        fprintf('UnKnownReduction /n')
    end
    
    if nargout > 1
        % computing the vectors the adjoint Born approximation acts on 
        z = -2*(husimi_mat')*(residual.*phi_u);
        
        % computing the rhs for the adjoint system
        rhs_adj = omega^2*z;
        
        % solving the adjoint system
        W_adj = (H1')\rhs_adj;
        
        % computing the gradient
        grad = real(sum(conj(U).*W_adj, 2))*h^2;
        
        % reshaping and extrating the gradient
        grad = reshape(grad, nx,ny);
        grad = grad(npml+1:npml+nxi, npml+1:npml+nyi);
        
        % returnig the gradient
        if reduction == "sum"
            dmis = grad(:)/(Ntheta_r*Ntheta_o);
        elseif reduction == "mean"
            dmis = grad(:)/(Ntheta_s*Ntheta_i*Ntheta_r*Ntheta_o);
        else
            fprintf('UnknownReduction')
        end
    
    end 
    
    
end
    