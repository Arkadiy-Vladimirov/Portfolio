function phasePortrait(ode,tspan,init_pnts,ode_name,varargin)
%PHASEPORTRAIT function draws phase portrait of ode given
    opts = odeset('Refine',10,'RelTol',1e-3,'AbsTol',1e-6); 
    arr_rar = 10; %rareness of velocity arrows;
    arr_sc = 1; %velocity arrows scale;
    
    figure('Name', 'phase potrait');
    hold on;
    for x0 = init_pnts
        [~,sol] = ode45(ode,tspan,x0,opts);
        x = sol(:,1); y = sol(:,2);
        if (abs(nargin) == 5)
            lyap_f = varargin{1};
            pl = colored_plot(x,y,lyap_f);
        else
            pl = plot(x,y,'Color','#0072BD');
        end
        
        sol = sol(1:arr_rar:end,:);
        w = ode(0,sol'); u = w(1,:); v = w(2,:);
        qu = quiver(sol(:,1),sol(:,2),u',v',arr_sc,'Color','#D95319');
    end
    ax = gca; 
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex'; ax.YLabel.Interpreter = 'latex';
    ax.XLabel.String = '$x$'; ax.YLabel.String = '$y$';
    ax.YLabel.Rotation = 0;
    ax.XLim = [-5 5]; ax.YLim = [-5 5];
    title(ode_name);
    
    if (abs(nargin) == 5)
        co = fcontour(lyap_f);
        legend([pl,qu,co],'phase curves','velocity',...
                            'lyapunov function levels');
    else
        legend([pl,qu],'phase curves','velocity');
    end
    
end

function plot_obj = colored_plot(x,y,l_f)
    N = length(x);
    col1 = [0.4940 0.1840 0.5560]; %lyap -> inf
    col2 = [0.3010 0.7450 0.9330]; %lyap -> 0
    alpha = @(x,y) -1./(l_f(x,y)+1) + 1;
    col_func = @(x,y) col1.*alpha(x,y) + col2.*(1-alpha(x,y));
    for cnt = 1:N-1
        col = col_func(x(cnt),y(cnt));
        x_seg = [x(cnt) x(cnt+1)]; y_seg = [y(cnt) y(cnt+1)]; 
        plot(x_seg,y_seg,'Color', col);
    end
    col = col_func(x(N),y(N));
    plot_obj = plot(x(N),y(N),'Color', col);
end