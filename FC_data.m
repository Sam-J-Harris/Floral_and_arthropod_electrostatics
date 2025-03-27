function [Vdt, VEdt, maxV, Evdt] = FC_data(zdt,V,VEswit,meswit,xspc,yspc)
% = applies e-potential V function handle to mesh and outputs e-potential and e-field data
Vdt = V(zdt); % compute e-potential
% E-field computation (if required)
if VEswit==0 % e-potential
    VEdt = Vdt; % e-pot dataset
    maxV = max(max(VEdt)); % max value of e-pot
    Evdt = 0; % e-field vectors (zero as not calculated here)
else % e-field magnitude
    if meswit ==0 % x and y
        [Ex, Ey] = gradient(-Vdt,xspc,yspc); % e-field vector 
    else % r and theta
        r = sqrt(real(zdt).^2+imag(zdt).^2); 
        [Ex, rmEy] = gradient(-Vdt,xspc,yspc); % e-field vector - x component
        Ey = rmEy./r; % e-field vector - y component
    end
    Emag = sqrt(Ex.^2+Ey.^2); % e-field magnitude
    VEdt = Emag; % e-mag dataset
    maxV = max(max(VEdt)); % max value of e-mag
    Evdt = 10.*[Ex, Ey, 0*ones(size(Ex))]; % e-field vector
end
end