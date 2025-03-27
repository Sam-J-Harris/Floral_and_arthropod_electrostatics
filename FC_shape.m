function z = FC_shape(cswit,Pts,npet,ang,A,B,C)
% = outputs the x,y coordinates of points on the flower boundary in the form z = x+iy.
scl = 1; rot = scl.*exp(1i*ang); % scaling (=1), rotation of the flower.
if cswit ==1 % pointy flower
    tht = linspace(-pi,pi,Pts); zta=exp(1i*tht); % azimuthal angle
    z = rot.*(A.*zta+B.*zta.^(-(npet-1))).'; % z coordinates
else    % smooth flower
    t = linspace(-1,1,Pts);
    z = (rot.*exp(1i*pi*t).*(1+C.*cos(npet*pi*t))).'; % z coordinates
end
z = z./max(abs(z)); % ensure a petal length of 1
end