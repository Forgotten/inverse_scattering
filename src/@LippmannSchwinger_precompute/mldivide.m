function x = mldivide(LS,b)
% MLDIVIDE  Overload the backslash operator for Ham class
%    X = H/B returns x that solevd the linear system
%
%  Copyright (c) 2019-2020 Leonardo Zepeda-Núñez
%  This file is distributed under the terms of the MIT License.

% x = gmres(@(u) LS*u, b, 20, 1e-6, 100);
[x, exitflag] = bicgstab(@(u) LS*u, b, 1e-8, 5000);

if exitflag ~= 0
    fprintf('BiCGStab failed to converge: %d /n', exitflag)
end

end
