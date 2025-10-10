%% backward_euler_step
t = 0;
h = 0.01;
% rate_func01
X0 = 1;

% rate_func01
X0 = 1;
forward_euler_step(@rate_func01, t, X0, h)
solution01(t+h)

backward_euler_step(@rate_func01, t, X0, h)


clear;
clc;
clf;
close all;

%% forward_euler_step
t = 0;
h = 0.3;

% rate_func01
X0 = 1;
forward_euler_step(@rate_func01, t, X0, h)
solution01(t+h)

% rate_func02
X0 = [1; 0];
forward_euler_step(@rate_func02, t, X0, h)
solution02(t+h)

%% forward_euler_fixed_step_integration
t_span = [0, 6];
h_list = [0.5, 0.4, 0.25];
colors = ['b', 'k', 'm', 'r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
    plot(t_list, X_list, colors(j));
    hold on;
    
    
end
ylim([-2,2])
xlabel('Time');
ylabel('Solution');
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(i) = solution01(t_list(i));
end
plot(t_list, solution_X_list, colors(4))
title('Forward Euler H Step Comparison');
legend("Euler's Method H = 0.5", "Euler's Method H = 0.4", "Euler's Method H = 0.25",'Closed-Form Solution')
grid on;



%% forward_euler_fixed_step_integration
t_span = [0, 10];

% rate_func01
X0 = 1;
[t_list, X_list, h_avg, num_evals] = fixed_step_integration(@rate_func01, @backward_euler_step, t_span, X0, h);
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(i) = solution01(t_list(i));
end

figure();
hold on;

plot(t_list, X_list, 'b.');
plot(t_list, solution_X_list)
xlabel('Time');
ylabel('Solution');

%% Defined functions
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end

function X = solution02(t)
    X = [cos(t);sin(t)];
end