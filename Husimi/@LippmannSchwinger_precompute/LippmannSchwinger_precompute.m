classdef LippmannSchwinger_precompute
    % TO DO: make this work in a rectangle! 
    % this class encodes the fast application of 
    % M = I + omega^2 G * spddiagm(nu)
    % all the properties are public
    
    % This is a version with precomputed Nystrom discretization
    
    properties (SetAccess = public)
        % array that contains the 
        GFFT
        % array nxm with the compactly supported heterogeneity
        nu
        % number of points in the extended domain
        ne
        me
        % number of points in the physical domain
        n
        m
        % frequency
        omega
        % domain
        x
        y
    end
     
    methods
        % constructor
        function LS = LippmannSchwinger_precompute(x,y,omega,nu)
            % x is a vector
            % y is a vector 
            % omega is a real
            % nu is either a function handle or a vector
            
            % by default we will use Greengard Vico
            
            % Domain parameters
            LS.n = length(x);
            LS.m = length(y);
            
            % two times extension for aperiodic convolution
            LS.ne = 2*LS.n;
            LS.me = 2*LS.m;
            LS.omega = omega;
            
            LS.x = x;
            LS.y = y;
            
            a = (x(end)-x(1));
            b = (y(end)-y(1));
            
            L  = ceil(sqrt(a^2+b^2)*2)/2; % oscillation in the truncated kernel
            
            % times of extension for initial step
            gamma = 4;
            % we need the computational domain in the first step at least 4 times bigger
            Lpx = gamma*a ;  
            Lpy = gamma*b ;
            
            % Frequency contents
            % n needs to be even
            kx = (-(gamma/2*LS.n):1:(gamma/2*LS.n-1));
            ky = (-(gamma/2*LS.m):1:(gamma/2*LS.m-1));

            KX = repmat((2*pi/Lpx)*kx',1,gamma*LS.m);
            KY = repmat((2*pi/Lpy)*ky,gamma*LS.n,1);

            % Truncated Green's function
            S = sqrt(KX.^2 + KY.^2);

            G2D = @(L,k,s)(1 + ...
                           (1i*pi/2*L*besselh(0,1,L*k)).*(s.*besselj(1,L*s)) - ...
                           (1i*pi/2*L*k*besselh(1,1,L*k)).*besselj(0,L*s)...
                           )./(s.^2 - k^2);

            GFFT_temp = fftshift(G2D(L, omega, S));

            
            % heterogeneity
%             LS.nu = reshape(nu(x',y),[],1);
            LS.nu = nu(:);
            % in some cases we (truncate the heterogeneity if it blows up)
            if any(isnan(LS.nu))
                fprintf("some values of the window may be Nan")
                LS.nu(isnan(LS.nu)) = 0.0;
            end
            
            % compute the fft of the Nystrom discretization
            % delta matrix
            T = zeros(gamma*LS.n,gamma*LS.m);
            T(1,1) = 1;
            % compute convolution matrix
            TFFT = fft2(T);
            TFFT = GFFT_temp.*TFFT;
            
            T = circshift(ifft2(TFFT),[LS.n,LS.m]);
            T = T(1:LS.ne,1:LS.me);
            
            % compute fft of convolution matrix
            LS.GFFT = fft2(fftshift(T));
        end
        
        % we need to use this one for the multiplication 
        %(input -> vector/matrix, output -> matrix)
        function Gu = apply_Green(LS, u)
            u = reshape(u,LS.n,LS.m);
            % Fourier Transform (zero padding)
            Gu = fft2(u,LS.ne,LS.me);
            % Component-wise multiplication
            Gu = LS.GFFT.*Gu;
            % Inverse Fourier Transform
            Gu = ifft2(Gu);
            % Extract useful information
            Gu = Gu(1:LS.n, 1:LS.m);
        end
        
%         % Plot
%         function plot(LS,u,BD,MD)
%             imagesc(LS.x,LS.y,u.'); pbaspect([1 1 1]);
%             ylabel('y'); xlabel('x'); set(gca,'YDir','normal'); colorbar;
%             if nargin >= 3
%                 % plot boundary
%                 hline_b = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color','k');
%                 rectangle('Position',[BD.xb-BD.rb BD.yb-BD.rb 2*BD.rb 2*BD.rb],...
%                     'Curvature',[1,1],'LineWidth',2);
%             end
%             if nargin >= 4
%                 % plot support
%                 rectangle('Position',[MD.xn-MD.rn MD.yn-MD.rn 2*MD.rn 2*MD.rn],...
%                     'Curvature',[1,1],'LineWidth',2,'LineStyle','--');
%                 hline_n = line(NaN,NaN,'LineWidth',2,'LineStyle','--','Color','k');
%                 
%                 % legend
%                 legend([hline_b,hline_n],'$\partial\Omega$','supp$(n-1)$','Interpreter','latex');
%             end
%         end
        
    end
end