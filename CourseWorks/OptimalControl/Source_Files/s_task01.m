%linear (processing speed?) problem
clear; clc;
load configs/config1.mat
t0 = 0;

sol = lpsp(A,B,f,t0,Pconf,r,pnts,params);

if (~isempty(sol))
disp(strcat("right transversality condition error ", num2str(sol.error)));
disp(strcat("optimal time ", num2str(sol.T)));
disp(strcat("(set reaching) phi interval ", num2str(sol.phi_int(1)),...
                                          " - ", num2str(sol.phi_int(2))));
exportgraphics(gcf,'examples/opt_traj_exmp.pdf','ContentType','vector');
plotControl(sol.t,sol.control);
exportgraphics(gcf,'examples/opt_cntrl_exmp.pdf','ContentType','vector');
plotConjugate(sol.conj_var);
exportgraphics(gcf,'examples/conj_exmp.pdf','ContentType','vector');

else
    disp("solution not found, try to increase maximum time")
    exportgraphics(gcf,'examples/traj_exmp.pdf','ContentType','vector');
end

%plotTrajectory(sol.solution.y);

%params.phi = sol.phi_int;
%sol = lpsp(A,B,f,t0,Pconf,r,pnts,params);
    
%% plotting functions
function plotControl(t,u)
    nameStr = 'Optimal control plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.Rotation = 0;
    ax.XLabel.String = '$t$';
    ax.YLabel.String = '$u_i$';
    hold on;
    
    plot(ax,t,u(1,:)); plot(ax,t,u(2,:)); 
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    legend('$u_1(t)$', '$u_2(t)$');
end

function plotTrajectory(x)
    nameStr = 'Optimal trajectory plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.Rotation = 0;
    ax.XLabel.String = '$x_1$';
    ax.YLabel.String = '$x_2$';
    hold on;
    
    plot(ax,x(1,:),x(2,:)); 
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    legend('$x(t)$');
end

function plotConjugate(psi)
    nameStr = 'Conjugate variable plot';
    fg = figure('Name', nameStr);
    ax = axes;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.Rotation = 0;
    ax.XLabel.String = '$\psi_1$';
    ax.YLabel.String = '$\psi_2$';
    hold on;
    
    plot(ax,psi(1,:),psi(2,:)); 
    lgd = legend(ax);
    lgd.Interpreter = 'latex';
    legend('$\psi(t)$');
end