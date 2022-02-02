function LSx = mtimes(LS,x)
% MTIMES  Overload multiplication operator for LippmannSchwinger class
%    LSx = LS*x returns a vector 
%
%    See also LippmannSchwinger

%  Copyright (c) 2019-2020 Leonardo Zepeda-Núñez

% we need to assert that x is a vector and that is has the correct size
assert(size(x,1) == LS.n*LS.m )
% only works for vectors so far
assert(size(x,2) == 1 )

% by construction this will give back a matrix
B = apply_Green(LS, x);

LSx = -x + LS.omega^2*LS.nu.*B(:);

end