function Ve = FC_exa(Z,tileps,efparam,Zpmat,Ajmat)
% = function handle for the e-potential V exact solution
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); %determines if point is in polygon
alpha = (tileps-1)/(tileps+1); beta = 1+alpha;

% e-field exact soln
V1e = @(z) -efparam.*beta.*abs(z).*cos(angle(z));
V2e = @(z) -efparam.*(abs(z).*cos(angle(z))+alpha.*cos(angle(z))./abs(z));

% arthropods exact soln
for j=1:size(Zpmat,2)
    zpt = Zpmat{j}; Aj = Ajmat{j};
    inout = inpolygonc(zpt,Z);

    if inout
        V1e = @(z) V1e(z)+Aj.*tileps.*(log(abs(z-zpt))-alpha.*log(abs(z.*conj(zpt)-1)));
        V2e = @(z) V2e(z)+Aj.*(beta.*log(abs(z-zpt))-alpha.*log(abs(z)));
    else
        V1e = @(z) V1e(z)+Aj.*(beta.*log(abs(z - zpt))-alpha.*log(abs(zpt)));
        V2e = @(z) V2e(z)+Aj.*(log(abs(z-zpt))+alpha.*log(abs(1./(z.*conj(zpt))-1)));
    end
end

% Combined exact soln
Ve = @(z) V1e(z).*inpolygonc(z,Z) + V2e(z).*(~inpolygonc(z,Z));
end