function f = polygon_bump(X,Y,cx,cy,r,p)

if round(p)~=p || (p <3 && p >0) || (p<0 && p~=-5 && p~=-6)
    error('p has to be zero, or positive integer not smaller than 3, or -5, -6.');
end

if p == 0                                      % circular support
    temp = 1 - ((X-cx).^2+(Y-cy).^2)/r^2;
    f = zeros(size(X));
    f(temp>0) = exp(1-1./temp(temp>0));
elseif p>0                                     % regular polygon support
    aa = atan2(Y-cy,X-cx)+pi/p;
    aa = mod(aa,2*pi/p)-pi/p;
    rr = r./(2*cos(aa));
    
    temp = 1 - ((X-cx).^2+(Y-cy).^2)./(rr.^2);
    f = zeros(size(X));
    f(temp>0) = exp(1-1./temp(temp>0));
elseif p<0                                     % star polygon
    p = abs(p);
    aa = atan2(Y-cy,X-cx)+pi/p;
    aa = mod(aa,2*pi/p)-pi/p;
    aa = abs(aa);
    rr = r*cos(2*pi/p)./cos(2*pi/p-aa);
    
    temp = 1 - ((X-cx).^2+(Y-cy).^2)./(rr.^2);
    f = zeros(size(X));
    f(temp>0) = exp(1-1./temp(temp>0));
end

end