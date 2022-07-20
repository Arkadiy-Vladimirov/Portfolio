%task03 script

%% №1
    clear; clc;
    N = 20;
    S = -(pi)^2/12;
    a_n = @(n)((-1).^n)./(n.^2);
    psi_n = @(n) 1./(n.^2);
    
    A1_N = a_n(1:N);
    S1_N = cumsum(A1_N);
    Y = abs(S1_N - S);
    PSI = psi_n(1:N);
    
    nameStr = 'series convergence plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = '$n$';
    hold on;
    plot(ax,Y); plot(ax,PSI);
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    legend('$|S_n - S|$','$\psi_n$');
    
%% №2
    clear; clc;
    N = 5; %number of attempts
    xint = [-3*pi/2, pi/2];
    lft_sd = @(x) cos(x);
    rgt_sd = @(x) x./pi;
    fun = @(x) lft_sd(x) - rgt_sd(x); 
    
    nameStr = 'roots of equation plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    fplot(ax, lft_sd, xint); hold on;
    fplot(ax, rgt_sd, xint);
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    labels = {'$\cos(x)$','$\frac{x}{\pi}$'};
    legend(labels);
    
    [x0_arr,~] = ginput; N = length(x0_arr);
    
    labels = cell(1,N+2);
    labels(1) = {'$\cos(x)$'}; labels(2) = {'$\frac{x}{\pi}$'};
    for cnt = 1:N
        x0 = x0_arr(cnt);
        [root,fval] = fzero(fun,x0); root
        misalignment = abs(0 - fval)
        pl = plot(root,lft_sd(root),'o');
        labels(cnt+2) = {strcat('root (',int2str(cnt),')')};
    end
    legend(labels);
    
    
%% №3
clear; clc;
    a = 4;
    xint = [-a,a];
    
    nameStr = 'roots from initial approximation plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    fplot(ax, @fNo3, xint,'MeshDensity',2);
    %X = linspace(-a,a,1000); Y = f(X); plot(X,Y);
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    label = {'$x\sin(\frac{1}{x})$ roots'};
    legend(label);

%% №4 %fix axes, slow animation, shift circle
clear; clc;
    R = 3; %constant(same as in function)
    alpha = 1.05; %reduce in velocity after each bounce
    tspan = [0 100]; %visible time (to detect if the ball stopped)
    pos0 = [0 2]; vel0 = [1 1]; %initial position and velocity
    rCircEvntFcn = @(t,x) CircEvntFcn(t,x,R); %func of circle boundary
    
    x0  = [pos0(1); vel0(1); pos0(2); vel0(2)];
    options = odeset('Events',rCircEvntFcn,'InitialStep', 1e-10);
    sol = ode45(@odefun,tspan,x0,options);
    t = sol.x;
    x = sol.y(1,:);
    y = sol.y(3,:);
    hasntStopped = ~isempty(sol.ie);
    
    nameStr = 'bouncing ball animation';
    figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    
    phi = linspace(0,2*pi,100);
    X = cos(phi)*R; Y = sin(phi)*R;
    plot(X,Y);
    hold on;
    
    while hasntStopped
        x0 = getInitPoint(sol,alpha,R);
        sol = ode45(@odefun,tspan,x0,options);
        x = cat(2,x,sol.y(1,:));
        y = cat(2,y,sol.y(3,:));
        t = cat(2,t,sol.x);
        hasntStopped = ~isempty(sol.ie);
    end
    
    comet3(ax,x,y,t,0.025);
    
%% №5 %analitic solution does not always exist
clear; clc;
    %a = 2/3; b = 4/3; g = 1; d = 1; x0 = 1.5; y0 = 1.5; tend = 20;
    a = 1; b = 0.01; g = 1; d = 0.02; x0 = 20; y0 = 20; tend = 12;
    %a = 0.03; b = 0; g = 1; d = 0.01; x0 = 20; y0 = 20; tend = 12; % b = 0!!!
 
    tspan = [0 tend];
    ode = @(t,y) L_V_ode(t,y,a,b,g,d);
    frst_int = @(x,y) d*x - g*log(x) + b*y - a*log(y);
    f_impl = @(x,y) frst_int(x,y) - frst_int(x0,y0);
    %analytical solution for integral curve (only if b = 0!!!)
    x_t = @(t) x0*exp(a*t); t_t = @(t) t;
    y_t = @(t) y0*exp(-g*t + (d*x0/a)*(exp(a*t)-1));
    
    opts = odeset('Refine',10,'RelTol',1e-3,'AbsTol',1e-6);
    [t,sol] = ode45(ode,tspan,[x0;y0],opts);
    x = sol(:,1); y = sol(:,2);
    
    fg1 = figure('Name', 'Lotka-Volterra equations');
    plot(x,y); hold on;
    fimplicit(f_impl);
    title('phase curve');
    xlabel('Prey population'); ylabel('Predator population'); 
    legend('numerical solution','analytical solution');
    
    fg2 = figure('Name', 'Lotka-Volterra equations');
    plot3(x,y,t); hold on;
    title('integral curve'); 
    xlabel('Prey population');ylabel('Predator population');zlabel('time');
    if (b == 0)
        fplot3(x_t,y_t,t_t);
        legend('numerical solution','analytical solution');
    else
        legend('numerical solution');
    end
    
%qualitative change in system behavior takes place if
% a/b/c/d = 0 (no cycles) or (x0,y0) = (g/d, a/b)

%% №6
clear; clc;
    N = 15;
    cnt_fnc = {@odeNode, @odeDicrNode, @odeSeddle, @odeFocus, @odeCentre;...
               'Node', 'Dicritic Node', 'Seddle', 'Focus', 'Centre'};
    tspan = [0 10];
    
    phi = linspace(0,2*pi,N); r = 10;
    init_pnts(1,:) = r.*cos(phi); init_pnts(2,:) = r.*sin(phi);
    
    ode = @odeNode; ode_name = 'Node';
    phasePortrait(ode,tspan,init_pnts,ode_name);
    ode = @odeDicrNode; ode_name = 'Dicritic Node';
    phasePortrait(ode,tspan,init_pnts,ode_name);
    ode = @odeFocus; ode_name = 'Focus';
    phasePortrait(ode,tspan,init_pnts,ode_name);
    
    phi = linspace(0,2*pi,2*N); r = 5;
    init_pntsS(1,:) = r.*cos(phi); init_pntsS(2,:) = r.*sin(phi);
    ode = @odeSeddle; ode_name = 'Seddle';
    phasePortrait(ode,tspan,init_pntsS,ode_name);
    
    init_pnts(1,:) = linspace(1,10,N); init_pnts(2,:) = zeros(1,N);
    ode = @odeCentre; ode_name = 'Centre';
    phasePortrait(ode,tspan,init_pnts,ode_name);
   
    
%% №7
clc; clear;
    N = 20;
    phi = linspace(0,2*pi,N); r = 5;
    init_pnts(1,:) = r*cos(phi); init_pnts(2,:) = r*sin(phi);
    ode = @ode1;
    ode_name = 'ODE 1';
    tspan = [0 10];
    f_lyap = @(x,y) x.^2 + y.^2;
    phasePortrait(ode,tspan,init_pnts,ode_name,f_lyap);
    
    phi = linspace(0,2*pi,N); r = 5;
    init_pnts(1,:) = r*cos(phi); init_pnts(2,:) = r*sin(phi);
    ode = @ode2;
    ode_name = 'ODE 2';
    tspan = [0 2];
    f_lyap = @(x,y) x.^2 + y.^2;
    phasePortrait(ode,tspan,init_pnts,ode_name,f_lyap);
    
%% №8 %resolve
clear; clc;
    N = 15;
    anSol = @(x) x.*0;
    xmesh  = linspace(0,pi/2,N);
    solinit = bvpinit(xmesh, @guess);
    sol = bvp4c(@bvpfun,@bcfun,solinit);
    
    figure('Name', 'bvp solution');
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex'; ax.YLabel.Interpreter = 'latex';
    ax.XLabel.String = '$x$'; ax.YLabel.String = '$y$';
    ax.YLabel.Rotation = 0;
    ax.XLim = [0 pi/2]; ax.YLim = [-3 3];
    hold on;
    fplot(ax,anSol); 
    plot(ax, sol.x,sol.y(2,:),'LineStyle','none','Marker','o');
    legend('analytical solution', 'numeric solution');
    
    L2_norm_difference = sqrt(trapz(sol.x,sol.y(2,:).^2))
    C_norm_difference = max(abs(sol.y(2,:)))
    
%% №9
clear; clc;
    %2arg function
    f = @f1; grad = @gradf1;
    x0 = [3; 5];
    sol = fmingd(f,x0,grad);
    
    figure('Name', 'gradient descent method optimization');
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex'; ax.YLabel.Interpreter = 'latex';
    ax.XLabel.String = '$x$'; ax.YLabel.String = '$y$';
    ax.YLabel.Rotation = 0;
    ax.XLim = [-10 10]; ax.YLim = [-10 10];
    hold on;
    
    fcont = @(x,y) f1([x y]);
    fcontour(ax,fcont,'LevelList',sol.stps_val);
    plot(ax, sol.stps(1,:),sol.stps(2,:),'Marker','.','MarkerSize',8);
    plot(ax, sol.xmin(1),sol.xmin(2),'Marker','o');
    legend('function levels', 'iterative steps','local minimum');
    
    %1 arg function
    f = @(x) x.^2; grad = @(x) 2*x;
    x0  = 3;
    my_solution = fmingd(f,x0,grad);
    [matlab_xmin matlab_fmin] = fminbnd(f,-10,10);
    my_solution.xmin
    matlab_xmin
    my_solution.fmin
    matlab_fmin
    
%% №10
clear; clc;
    fg1 = figure;
    %fg2 = figure;  %bug
    stp = 0.001; inpLims = [-10 10]; outLims = [-100 100];
    %info = plotFT(fg1, @testf2, @testft2, stp, inpLims, outLims)
    info = plotFT(fg1, @func1, @ftfunc1, stp, inpLims, outLims)
    fg2 = figure;
    info = plotFT(fg2, @func2, [], stp, inpLims, outLims);
    info = plotFT(fg2, @func2, @ftfunc2, stp, inpLims, outLims);
    %fg3 = figure;
    %info = plotFT(fg3, @func3, [], stp, inpLims, outLims);
    %fg4 = figure;
    %info = plotFT(fg4, @func4, [], stp, inpLims, outLims);
%% №11
clear; clc;
%in LaTeX.




%% functions block
%for №3
function root = fNo3(x0)
    N = length(x0);
    root  = zeros(N);
    for cnt = 1:N
        root(cnt) = fzero(@insinc,x0(cnt));
    end
    
    function y = insinc(x)
        if (x ~= 0)
            y = x.*sin(1./x);
        else
            y = 0;
        end
    end
end

%for №4
function [value,isterminal,direction] = CircEvntFcn(t,x,R)
    value = x(1)^2 + x(3)^2 - R^2;
    isterminal = 1;
    direction = 11;
end
function dxdt = odefun(t,x)
    dxdt = zeros(4,1);
    dxdt(1) = x(2);
    dxdt(2) = 0;
    dxdt(3) = x(4);
    dxdt(4) = 0;
end
function x0 = getInitPoint(solution,alpha,R)
    r1 = solution.ye(1); r2 = solution.ye(3);
    u1 = solution.ye(2); u2 = solution.ye(4);
    k = 1/(alpha*R^2);
    w1 = k*(u1*(r2^2 - r1^2) - 2*r1*r2*u2);
    w2 = k*(-2*r1*r2*u1 - u2*(r2^2-r1^2));
    x0 = [r1; w1; r2; w2];
end

%for №5
function dydt = L_V_ode(t,y,alpha,beta,gamma,delta)
    x = y(1); y = y(2); dydt = zeros(2,1);
    dydt(1) = alpha*x - beta*x*y;
    dydt(2) = -gamma*y + delta*x*y;
end

%for №6
function dydt = odeNode(t,y)
    A = [-2 -1; -1 -4];
    dydt = A*y;
end
function dydt = odeDicrNode(t,y)
    A = [-1 0; 0 -1];
    dydt = A*y;
end
function dydt = odeSeddle(t,y)
    A = [3 -2; 4 -6];
    dydt = A*y;
end
function dydt = odeFocus(t,y)
    A = [1 3; -6 -5];
    dydt = A*y;
end
function dydt = odeCentre(t,y) 
    A = [-2 -5; 2 2];
    dydt = A*y;
end

%for №7
function dzdt = ode1(t,z)
    dzdt = []; x = z(1,:); y = z(2,:);
    dzdt(1,:) = y - x + x.*y;
    dzdt(2,:) = x - y - x.^2 - y.^3;
end
function dzdt = ode2(t,z)
    dzdt = []; x = z(1,:); y = z(2,:);
    dzdt(1,:) = x.^2 + 2*y.^3;
    dzdt(2,:) = x.*y.^2;
end

%for №8
function dydx = bvpfun(x,y)
    dydx = [y(2)
            -y(1)+1];
end
function res = bcfun(ya,yb)
    res = [ya(2)
           yb(2)];
end
function y = guess(x)
    y = [2*cos(2*x)
        sin(2*x)];
end

%for №9
function fun = f1(x)
    fun = sum(x.^2);
    %fun  = sin(x(1))*cos(x(2)); %another example
end
function grad = gradf1(x)
    grad = [2*x(1), 2*x(2)];
    %grad = [cos(x(1))*cos(x(2)), -sin(x(1))*sin(x(2))]; %another example
end

%for №10
function f = func1(t)
    f = t.*exp(-2*abs(t)).*sinh(t);
end
function f = func2(t)
    f = t./(2 + 2*t + t.^2);
end
function f = func3(t)
    f = exp(-t.^6).*atan(t.^2);
end
function f = func4(t)
    mask = (abs(t+2) <= 1);
    t = t.*mask;
    f = atan(t.^3);
end
function f = testf1(t)
    f = 1./(t.^2 + 1);
end
function f = testf2(t)
    f = exp(-abs(t));
end

function F = ftfunc1(l) 
    F = 35*l.*(5+l.^2) ./ ((1+l.^2).^2 .* (9+l.^2).^2);
end
function F = ftfunc2(l)
    maskGZ = (l>0); maskEZ = (l==0); maskLZ = (l<0);
    F = pi*(-1-i)*exp(l.*(-1+i)).*maskGZ + ...
        + -pi.*maskEZ + ...
        + pi*(-1+i)*exp(l.*(1+i)).*maskLZ;
end
function F = testft1(l)
    F = pi*exp(-abs(l));
end
function F = testft2(l)
    F = 2./(1+l.^2);
end