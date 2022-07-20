function optim_sol = lpsp(A,B,f,t0,Pconf,r,pnts,params)
%LPSP function solves linear processing speed problem
%   Detailed explanation goes here

%TO DO:
%implement universal plotting function
%LaTeX

    %set constants
    nica = params.N; %size of initial conditions sample;
    dt = params.dt;  %time accuracy
    t1 = params.maxT;%max time
    tgrid = t0:dt:t1;
    
    %set auxilary functions
    ellX0_conf = {r,r,0,[0;0]};
    rhoX0 = @(psi) rhoEll(psi,ellX0_conf);
    ellP_conf = {sqrt(Pconf.alpha*Pconf.c),sqrt(Pconf.beta*Pconf.c),...
                    0,[Pconf.a;Pconf.b]};
    rhoP = @(psi) rhoEll(psi,ellP_conf);
    
    %auxilary data for event function(normals,vectors of sides)
    k = convhull(pnts,'Simplify',true); k = k(1:end-1); lk = length(k);
    pnts = pnts(k,:);
    vects = circshift(pnts,-1,1) - pnts;
    normals = ones(lk,2);
    normals(:,1) = vects(:,2); normals(:,2) = - vects(:,1);
   
    options = odeset('Events',@(t,x)inX1(t,x,pnts,normals));
    
    %initial conditions sample
    phi_int = params.phi;
    phi = linspace(phi_int(1), phi_int(2), nica); % 1xN
    psi0arr  = [cos(phi) ; sin(phi)]; %2xN
    [~,x0arr] = rhoX0(psi0arr); %2xN
    
    %plot preparation
    nameStr = 'Optimum suspicious trajectories plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex'; ax.YLabel.Interpreter = 'latex';
    ax.YLabel.Rotation = 0;
    ax.XLabel.String = '$x_1$'; ax.YLabel.String = '$x_2$';
    hold on;
    
    %draw set X0,X1
    X = r*cos(phi); Y = r*sin(phi);
    plX0 = plot(X,Y,'Color','#0072BD'); hold on;
    draw_pnts = pnts;
    k = convhull(draw_pnts,'Simplify',true);
    draw_pnts = draw_pnts(k,:);
    X = draw_pnts(:,1); Y = draw_pnts(:,2);
    plX1 = plot(X,Y,'Color','#D95319');
   
    %----------------------------------------------------------------------
    %SOLVING PROBLEM
    phi_opt = []; opt_sol = [];
    minT = t1;
    for cnt1 = 1:nica
        psi0 = psi0arr(:,cnt1);
        x0 = x0arr(:,cnt1);
        [~,psit] = ode45(@conjproblem,tgrid,psi0); psit = psit'; %2xN
        ut = getUt(tgrid,psit); 
        sol = ode45(@(t,x) ode(t,x,tgrid,ut),[t0 t1], x0, options);
        
        if (~isempty(sol.xe))
        if (sol.xe(end) <= minT)
            minT = sol.xe(end);
            x_T = sol.ye(:,end);
            opt_sol = sol;
            opt_ut = ut;
            opt_psi = psit;
            psi_T = zeros(2,1);
                psi_T(1) = interp1(tgrid,opt_psi(1,:),minT);
                psi_T(2) = interp1(tgrid,opt_psi(2,:),minT);
            err = abs(dot(-psi_T,x_T) - rhoX1(-psi_T, pnts));
            phi_opt = [phi_opt,phi(cnt1)]; % save all initial phi's with time less than t1
        end
        end
        %draw all alternatives
        plsmp = plot(sol.y(1,:),sol.y(2,:),'Color','#77AC30');
    end
    %----------------------------------------------------------------------
    
    if (~isempty(opt_sol))
        plopt = plot(opt_sol.y(1,:),opt_sol.y(2,:),'LineWidth',1.5,...
            'Color','#A2142F'); %draw the best one
        lg = legend([plX0,plX1,plsmp,plopt], {'$\partial\mathcal{X}_0$',...
            '$\partial\mathcal{X}_1$','sample trajectories','optimal trajectory'});
        lg.Interpreter = 'latex';
        lg.Location = 'southeast';
    
    %forming output
        t_len = length(tgrid(tgrid < minT));
        optim_sol.solution = opt_sol;
        optim_sol.control = opt_ut(:,1:t_len);
        optim_sol.conj_var = opt_psi(:,1:t_len);
        optim_sol.T = minT;
        optim_sol.t = tgrid(1:t_len);
        optim_sol.error = err;
        optim_sol.phi_int = [min(phi_opt) max(phi_opt)];
    else
        optim_sol = [];
        lg = legend([plX0,plX1,plsmp], {'$\partial\mathcal{X}_0$',...
            '$\partial\mathcal{X}_1$','sample trajectories'});
        lg.Interpreter = 'latex';
        lg.Location = 'southeast';
    end
    
%__________________________________________________________________________
%LOCAL FUNCTIONS SECTION

%conjugate variable ODE
function dpsidt = conjproblem(t,psi)
    dpsidt = - (A(t))' * psi;
end

%direct variable ODE
function dxdt = ode(t,x,tgrid,ut)
       u = zeros(2,1);
       u(1) = interp1(tgrid,ut(1,:),t);
       u(2) = interp1(tgrid,ut(2,:),t);
       dxdt = A(t)*x + B(t)*u + f(t);
end

%u(t) control function
function ut = getUt(tgrid,psi)
    N = length(tgrid);
    ut = zeros(2,N);
    for cnt2 = 1:N
        Bt = B(tgrid(cnt2));
        [~,suppnt] = rhoP(transpose(Bt)*psi(:,cnt2));
        ut(:,cnt2) = suppnt;
    end
end

%X1 set indicator function
function [value,isterminal,direction] = inX1(t,x,points,normals)
    rv = [x(1) x(2)]; N = length(points);
    value = 0;
    for i = 1:N
        if (dot(rv-points(i,:),normals(i,:)) > 0)
            value = 1; %out of X1
        end
    end
    isterminal = 1;
    direction = 0;
end

%X1 support function
function val = rhoX1(psi, pnts)
    n = length(pnts);
    vals = zeros(n,1);
    for cnt = 1:n
        vals(cnt) = dot(psi,pnts(cnt,:));
    end
    val = max(vals);
end

end


