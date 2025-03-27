function FC_plot(Vdata,VEdata,maxVE,Evec,errdata,meanerr,VEswit,COswit,errswit,z,efparam,Zpmat,Ajmat,xdt,ydt,fignum)
% = produces contour plots of e-potential/e-field magnitude.
%
% Code:
% Graph parameters
titlelist = {'AAA-LS','Exact'}; VEstr = " E-potential "; % potential graph titles
LW = 'LineWidth'; lwsz = 1.8; % width of plot lines

X = xdt; Y = ydt; 
xmin = min(min(X)); xmax = max(max(X)); ymin = min(min(Y)); ymax = max(max(Y));
axlst = [xmin xmax ymin ymax]; % axis data

xp = real(z); yp= imag(z); zp = (1e8).*ones(size(xp)); % flower shape data
zdt = 1e8*ones(size(X)); % place Evec arrows above surface plot
rg = 1:10:size(X,1); % range of Evec arrows plotted

Vend = size(VEdata,2); ECend = size(errdata,2); 
subrow = 1; subclm = Vend + +ceil(errswit/2); % total number of plots (1 regularly, 3 for error data) 

% Surface plot (and contour lines) of e-potential/e-magnitude
for k=1:Vend 
    figure(fignum)
    subplot(subrow,subclm,k)
    
    Z = VEdata{k}; % e-pot/e-mag data

    % Control contour limits (change as you like)           
    if VEswit==1 % for e-magnitude
        sftol = 5; Z(Z>sftol) = sftol; % upper limit of contours 
    elseif efparam ==0 && size(Ajmat,2)>0 % arthropods only (w/ approx charge=1)
        sftol = 0.1; Z(Z>sftol) = sftol;
    end

    % Surface plot
    surf(X,Y,(1-12)*ones(size(Z)),Z); hold on, % surface plot of e-pot/e-mag
    view(2); shading interp, 

    % Contour plot (if applicable)
    if COswit==1
        V = Vdata{k};
        contour(X,Y,V,50,'w'); % contours placed above surf plot
        %contourf(X,Y,Z+max(max(Z)).*ones(size(Z)),20); % filled contour plot contours placed above surf plot
    end

    % Flower shape and arthropods (uncomment if desired)
    % plot3(xp,yp,zp,'r',LW, lwsz); % plot flower shape
    % for j=1:size(Ajmat,2) % plot bees/spiders
    %     zpt = Zpmat{j}; Aj = Ajmat{j};
    %     if Aj~=0
    %         plot3(real(zpt),imag(zpt),zp,'*','Color','k','MarkerSize', 12, LW, lwsz);
    %     end
    % end

    % E-field vectors (if applicable)
    if VEswit==1
        %colormap(hot), 
        colormap jet
        VEstr = " E-magnitude "; % change title string to say E-mag
        qdt = Evec{k}; q1d = qdt(:,1); q2d = qdt(:,2); q3d = qdt(:,3); % e-field data
        %quiver3(X(rg),Y(rg),zdt(rg),q1d(rg),q2d(rg),q3d(rg),'w',LW,1.25); % plot quivers (comment on/off as desired)
    else
        colormap jet
    end
    
    % Graph properties
    colorbar, 
    title({titlelist{k} + VEstr + "Solution"}, 'FontSize', 16)
    subtitle({"Max value = " + num2str(maxVE{k})}, 'FontSize',14)
    daspect([1 1 1]), axis(axlst); axis square
end

% Error plots (only if errswit~=0)
if errswit~=0  
    for k=1:ECend
        figure(fignum)
        subplot(subrow,subclm,subclm)

        Z2 = errdata{k};

        surf(X,Y,Z2);
        view(2); shading interp, hold on, 
        plot3(xp,yp,zp,'r',LW, lwsz); % flower shape

        colorbar,
        title({titlelist{k} + " Error"}, 'FontSize',16)
        subtitle({"Max error = " + num2str(meanerr{k})}, 'FontSize',14)
        daspect([1 1 1]), axis(axlst);
    end
end

% Main title
if efparam~=0 && size(Ajmat,2)>0
    bigtitle = 'Uniform e-field and arthropods';
elseif efparam~=0
    bigtitle = 'Uniform e-field only';
else
    bigtitle = 'Arthropods only';
end
sgtitle(bigtitle, 'FontSize', 24, 'FontWeight','bold')
set(gcf, 'WindowState', 'maximized'); % final figure is maximised
end