function Va = FC_aaa(Z,tileps,efparam,Zpmat,Ajmat)
% = function handle for the e-potential V calculates by the two-domain AAA-least squares algorithm.
% Va computes both interior and exterior e-potential.
%
% Code:
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); %determines if point is in polygon

% Normal vector
dZ = circshift(Z(1:end-1),1) - circshift(Z(1:end-1),-1); v = 1i.*(dZ)./abs(dZ); 
nvec = v./abs(v); nvec(end+1) = nvec(1); % outward pointing normal vector for each point on the flower boundary

% Far field condition into Dirichlet (dBC) and Neumann (nBC) BCs - background (b.)e-field and arthropods (a.pods)
hd = @(z) -efparam.*real(z); % Dirichlet BC - b.e-field term
hn = @(z) -efparam.*tileps.*real(nvec); % Neumann BC - b.e-field term

V1term = @(z) 0.*z; % interior potential V1 - b.e-field term
V2term = @(z) -efparam.*real(z); % exterior potential V2 - b.e-field term

for j=1:size(Zpmat,2)
    zpt = Zpmat{j}; Aj = Ajmat{j};
    inout = inpolygonc(zpt,Z);  % determines if pt charge in/out of flower

    hd = @(z) hd(z)+Aj.*(1-tileps.*inout).*log(abs(z-zpt)); % Dirichlet BC - a.pods term
    hn = @(z) hn(z)+Aj.*tileps.*(1-inout).*real(nvec./(z-zpt)); % Neumann BC - a.pods term

    V1term = @(z) V1term(z) + Aj.*tileps.*inout.*log(abs(z-zpt)); % interior potential V1 - a.pods term
    V2term = @(z) V2term(z) + Aj.*log(abs(z-zpt)); % exterior potential V2 - a.pods term
end

Hd = hd(Z); Hn = hn(Z); bigH = [Hd; Hn]; % dBC and nBC applied to boundary data

% Global AAA poles %%%%%%%% 
rtoly = 1e-8; ztoly = 1e-2; % residue and distance tolerances
[~,poltot1,restot1] = aaa(Hd,Z,'cleanup',1,'toly',rtoly); % poles for Dirichlet BC
[~,poltot2,restot2] = aaa(Hn,Z,'cleanup',1,'toly',rtoly); % poles for Neumann BC

poltot = [poltot1; poltot2].'; restot = [restot1; restot2].'; % all poles and resolutions

% Pole control - manual elimination of Froissart doublets  %%%%%%%%
poltot = poltot(abs(restot)>rtoly); % poles of small residue < rtoly eliminated
D = min(abs(Z-poltot),[],1); D(D<ztoly*max(abs(Z))) = 0; D(D~=0) =1; % finds if min dist. b/ween poles and bdry is below tolerance
poltot=poltot.*D; poltot(abs(poltot)==0)=[]; % poles close to the boundary < ztoly eliminated

polint = (poltot(~inpolygonc(poltot,Z))); % exterior poles (used in interior problem)
polext = (poltot(inpolygonc(poltot,Z))); % interior poles (used in exterior problem)

% % Plot of poles %%%%%%%%
% figure(2) % uncomment if you want to see figure of the poles
% zplot = Z(1:end-1); xp = real(zplot); yp = imag(zplot);
% plot(xp,yp,'k'), hold on, plot(polint,'.r'); plot(polext,'.b'), hold off, 
% axis square, daspect([1 1 1]); axis([-10 10 -10 10]);

% Interior and exterior poles with min distance to bdry
if size(polint,1)==0 && size(polint,2)==0 % interior poles
    polint=0; dint=0; % if no poles produced, set these variables to zero
else
    dint = min(abs((Z)-polint),[],1); % d=distance between pole and closest point on polygon
end

if size(polext,1)==0 && size(polext,2)==0 % exterior poles
    polext=0; dext=0; % if no poles produced, set these variables to zero
else
    dext = min(abs((Z)-polext),[],1); % d=distance between pole and closest point on polygon
end

% Least-Squares %%%%%%%% 
N = 20; N2int = size(polint,2); % truncate power series at N
P1d = Z.^(1:N); % interior Dirichlet Runge bases
P2d = (1./Z).^(1:N); % exterior Dirichlet Runge bases

P1n = (1:N).*Z.^(0:N-1).*nvec; % interior Neumann Runge bases
P2n = -(1:N).*Z.^(-(2:N+1)).*nvec; % exterior Neumann Runge bases

Q1d = dint./((Z-polint)); % interior Dirichlet Neumann bases
Q2d = dext./((Z-polext)); % exterior Dirichlet Neumann bases

Q1n = -dint./((Z-polint).^2).*nvec; % interior Dirichlet Neumann bases
Q2n = -dext./((Z-polext).^2).*nvec; % exterior Dirichlet Neumann bases
      
A1d = [real(P1d) real(Q1d) -imag(P1d) -imag(Q1d)];
A2d = [real(P2d) real(Q2d) -imag(P2d) -imag(Q2d)]; 

A1n = [real(P1n) real(Q1n) -imag(P1n) -imag(Q1n)];
A2n = tileps.*[real(P2n) real(Q2n) -imag(P2n) -imag(Q2n)];

bigA = [A1d A2d; A1n A2n]; % big matrix of basis vectors for both dBC and nBC

testc = bigA\bigH; % least squares procedure to find unknown coefficients in power series

c1 = testc(1:(2*N+2*N2int),:); c2 = -testc((2*N+2*N2int)+1:end,:);
c1 = reshape(c1,[],2)*[1;1i]; c2 = reshape(c2,[],2)*[1;1i]; % coefficients c1 for V1, c2 for V2

F1 = @(z) reshape([z(:).^(1:N) dint./(z(:)-polint)]*c1, size(z)); % analytic function F1(z)
F2 = @(z) reshape([(1./z(:)).^(1:N) dext./(z(:)-polext)]*c2, size(z)); % analytic function F2(z)

V1a = @(z) V1term(z) + real(F1(z)); % interior potential V1
V2a = @(z) V2term(z) + real(F2(z)); % exterior potential V2

Va = @(z) V1a(z).*inpolygonc(z,Z) + V2a(z).*(~inpolygonc(z,Z)); % full potential function handle
end