function [tileps,efparam,Zpmat,Ajmat] = FC_nondim(L,Einf,lambmat,xpmat,ypmat,z,eps0,eps1,eps2)
% = outputs non-dimensional parameters, each is defined below.
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); %determines if point is in polygon (ie inside the flower)

    tileps = eps2/eps1; % ratio of permittivities: eps2/eps1 = exterior/interior permittivity
    efparam = (Einf~=0); % e-field param: determines if there is a background e-field (b.e-field)
    Zpmat = {}; Qmat = {}; Ajmat = {}; % cell arrays for pt charge data
    for j = 1:size(lambmat,2)
        zp = xpmat{j} + 1i*ypmat{j}; Zpmat{j} = zp; % position of arthropod
        inout = inpolygonc(zp,z); % determines if arthropod inside/outside of flower
        Qmat{j} = -lambmat{j}./(4*pi*eps0.*(eps2*(1-inout)+eps1*inout)); % Q = lambda/2*pi*eps
        Ajmat{j} = Qmat{j}./(L*Einf*efparam+Qmat{1}*(1-efparam)); % no b.e-field: A_j = Q_j/Q_1, w/ b.e-field: A_j = Q_j/(L*E_inf) 
    end
end