function [xdt, ydt, zdt, xspc, yspc] = FC_mesh(M,N,Rmax,meswit)
% = outputs points in the mesh grid that AAA-LS will solve across. 
if meswit==0 % Cartesian grid
    nX = linspace(-Rmax,Rmax,M); nY = linspace(-Rmax,Rmax,N);
    [xdt,ydt] = meshgrid(nX,nY); zdt = xdt+1i*ydt; % creates xdt,ydt meshgrid and creates complex variable zdt = xdt+iydt
    xspc = nX(2)-nX(1); yspc = nY(2)-nY(1); % spacing between adjacent x and y points
else % Radial grid
    nR = linspace(0,Rmax,M); nT = linspace(0,2*pi,N); 
    [R,T] = meshgrid(nR,nT); xdt = R.*cos(T); ydt = R.*sin(T); zdt = xdt+1i*ydt; % creates radial R,T meshgrid, then creates xdt, ydt and zdt from that
    xspc = nR(2)-nR(1); yspc = nT(2)-nT(1); % spacing between adjacent x and y points
end
end