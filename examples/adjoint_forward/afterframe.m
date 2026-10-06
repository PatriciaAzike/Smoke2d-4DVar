%parms = read_vars();
%W = parms.W_eps;

title(sprintf('$\\hat{q}$ at t = %.2f', t), 'Interpreter', 'latex')
%title(sprintf('$\\hat{q}$ at time %.2f', t), 'Interpreter', 'latex')

pseudo_1d = false;

pseudo_1d_data = 'Uhat_data.mat';
plot_soln = false;
plot_surf = true;
plot_1d_data = false;


fprintf("qmin = %12.4e\n",qmin);
fprintf("qmax = %12.4e\n",qmax);


if PlotType ==1
    if plot_soln
        N = 500;
        if (t > 0)
            [xout,yout] = filament_soln(N,t);
        else
            th = linspace(0,2*pi,N+1);
            xout = 0.25*cos(th) + 0.5;
            yout = 0.25*sin(th) + 1;
        end
        hold on;
        plot(xout,yout,'k','linewidth',2);
        hold off;
    end



    % Colormap and axis
    colormap("parula"); %yrbcolormap
    clim([-0.5,4.5])
    cb = colorbar;  
    cv = linspace(qmin,qmax,11);
    cv(1) = [];
    %drawcontourlines(cv);

    % grids and borders
    showpatchborders(1:5);
    setpatchborderprops('linewidth',1)
    set(gca,'FontSize',16)

    % axes
    axis([0 2 0 2])
    %daspect([1 1 1]);
    view(2);
    %shg;

    if ~pseudo_1d
        hdl = add_gauges();
        set(hdl,'markersize',20);
        set(hdl,'color','k')
        ud = get(hdl,'Userdata');
        %set(ud{1},'string','');
        %set(ud{2},'string','');     
    end

   
    xlabel('x','FontSize',16);
    ylabel('y','FontSize',16);

    style_panel(cb, 0:4);
    shg;

    frame_dir = fullfile(fileparts(mfilename('fullpath')), 'amr_frames');
    if ~exist(frame_dir, 'dir')
        mkdir(frame_dir);
    end
    frame_file = fullfile(frame_dir, sprintf('amr_%04d_t%05.2f.png', Frame, t));
    print(gcf, frame_file, '-dpng', '-r300');
    % if exist('exportgraphics', 'file') == 2
    %     exportgraphics(gcf, frame_file, 'Resolution', 300);
    % else
    %     saveas(gcf, frame_file);
    % end

    if plot_surf
        figure(2)
        clf;
        [xcm,ycm] = meshgrid(xcenter,ycenter);
        surf(xcm,ycm,q','edgecolor','none');

        hold on;
        plot3(0.67, 1, 3, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
        plot3(0.8, 1.2, 3, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
        plot3(0.55, 1.2, 2, 'kd', 'MarkerSize', 15, 'MarkerFaceColor', 'k');

        set(gca, 'ZLim', [0 5]);
        camlight;
        title(sprintf('t = %.2f', t), 'Interpreter', 'latex')
        set(gca,'FontSize',16)

        figure(1);
    end




    %[tcurr_hdl,track_hdl] = track_gauges(t,'all');

    NoQuery = 0;
    prt = false;
    if (prt)
        hidepatchborders;
        figsize=[4,4];
        minlevel = 3;
        maxlevel = 5;
        mi = 2;
        mj = 2;
        mx = 8;
        maxres = mi*mx*2^maxlevel;   % 8*32 = 256
        dpi = maxres/4;
        plot_tikz_fig(Frame,figsize,'plot',dpi);
    end

elseif PlotType == 4 & plot_1d_data
    hclaw = findobj(gca, '-property', 'XData');

    %figure('Position',[100 100 1000 700]);
    axis([0,2,-0.1,4]);
    hold on;
    x_F = linspace(0,2,1024);
    
    %q_F = (abs(mod(x_F-t-0.5,2)) < 0.25);
    d   = mod(x_F - t - 0.5 + 1, 2) - 1;   % periodic difference in [-1, 1)
    q_F = (abs(d) < 0.25);
    h=plot(x_F,q_F,'r.');
    %plot(0.67,3,'rv','MarkerSize',20);
    h3=plot(0.67,3,'r*','MarkerSize',15,'LineWidth',2);
    h4=plot(0.8,3,'r*','MarkerSize',15,'LineWidth',2);
    h5=plot(0.55,2,'r*','MarkerSize',15,'LineWidth',2);
    %h6=plot(0.6,1,'r*','MarkerSize',20);

    % PLot data from 1D_advection code with representers in python
    X = load(pseudo_1d_data);

    tv = X.tv_hat;

    n = max(find(tv<=t));                 
    U_hat = X.Uhat(:,n); 
    
    %h1 = plot(X.xc_hat, U_hat, 'm-');

    xlabel('x','FontSize',16)
    ylabel('PM$_{2.5}$ concentration','FontSize',16,'Interpreter','latex')
    
    set(gca,'FontSize',16)

    h2 = getlegendinfo();
    h2 = h2(ishghandle(h2));

    
    % legend_handles = [h,h1,h3,h4,h5];
    legend_handles = [h,h3,h4,h5];
    legend_labels = {'$q_F$', 'obs1', 'obs2', 'obs3'};
    % legend_labels = {'$q_F$', '1D Notebook $\hat{q}$', 'obs1', 'obs2', 'obs3'};

    if ~isempty(h2)
        legend_handles = [h,h2(1),h3,h4,h5];
        legend_labels = {'$q_F$', '$\hat{q}$', 'obs1', 'obs2', 'obs3'};
    end

    legend(legend_handles, ...
        legend_labels, ...
        'Interpreter','latex', ...
        'FontSize',16, ...
        'Location','northeast')
    

    set(gca,'Box','on')

    frame_dir = fullfile(fileparts(mfilename('fullpath')), 'pseudo-1d');
    if ~exist(frame_dir, 'dir')
         mkdir(frame_dir);
    end
    frame_file = fullfile(frame_dir, sprintf('pseudo_%04d_t%05.2f.png', Frame, t));

    
    thin_scatter(hclaw);

    style_line_panel();
    print(gcf, frame_file, '-dpng', '-r300');

    hold off;
   
end

shg
