clear;
clc;
clf;
close all;

%% forward_euler_step
t = 0;
h = 0.01;

% rate_func01
X0 = 1;
forward_euler_step(@rate_func01, t, X0, h)
solution01(t+h)

% rate_func02
X0 = [1; 0];
forward_euler_step(@rate_func02, t, X0, h)
solution02(t+h)

%% forward_euler_fixed_step_integration
t_span = [0, 100];

% rate_func01
X0 = solution01(0);
[t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, t_span, X0, h);
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(i) = solution01(t_list(i));
end

figure();
hold on;

plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time');
ylabel('Solution');
title('Test01 Forward Euler Integration Results');
legend("Euler's Method", 'Closed-Form Solution')
grid on;

% rate_func02
X0 = solution02(0);
[t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func02, t_span, X0, h);
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(:, i) = solution02(t_list(i));
end

figure();
hold on;

plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time');
ylabel('Solution');
title('Test02 Forward Euler Integration Results');
grid on;

%Midpoint

% rate_func01
X0 = solution01(0);
[t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h);
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(i) = solution01(t_list(i));
end

figure();
hold on;

plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time');
ylabel('Solution');
title('Test01 Explicit Midpoint Integration Results');
legend("Midpoint Method", 'Closed-Form Solution')
grid on; 

% rate_func02
X0 = solution02(0);
[t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func02, t_span, X0, h);
solution_X_list = zeros(length(X0), length(t_list));
for i = 1:length(t_list)
    solution_X_list(:, i) = solution02(t_list(i));
end

figure();
hold on;

plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time');
ylabel('Solution');
title('Test02 Explicit Midpoint Integration Results');
grid on;

%% Local Truncation Error
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, 1, num_points);
X0 = solution01(t_ref);
x_actual_list = zeros(1, length(h_list));

% Euler's method
euler_x_list = zeros(1, length(h_list));
for i = 1:num_points
    [X1, num_evals] = forward_euler_step(@rate_func01, t_ref, X0, h_list(i));
    euler_x_list(i) = X1;
    x_actual_list(i) = solution01(t_ref + h_list(i));
end

% Midpoint method
midpoint_x_list = zeros(1, length(h_list));
for i = 1:num_points
    [X1, num_evals] = explicit_midpoint_step(@rate_func01, t_ref, X0, h_list(i));
    midpoint_x_list(i) = X1;
end

% Calculate global truncation error
euler_truncation_error = abs(x_actual_list - euler_x_list);
midpoint_truncation_error = abs(x_actual_list - midpoint_x_list);

figure();
loglog(h_list, euler_truncation_error)
hold on;
loglog(h_list, midpoint_truncation_error)
xlabel('Step Size (h)');
ylabel('Local Truncation Error');
title('Local Truncation Error vs Step Size');
legend('Euler Error', 'Midpoint Error');
grid on;

%% Global Truncation Error
X0 = solution01(0);
t_span = [0, t_ref];
x_actual = solution01(t_ref);

% Euler's method
euler_x_list = zeros(1, length(h_list));
for i = 1:num_points
    [t_list, X_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, t_span, X0, h_list(i));
    euler_x_list(i) = X_list(:, end);
end

% Midpoint method
midpoint_x_list = zeros(1, length(h_list));
for i = 1:num_points
    [t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_list(i));
    midpoint_x_list(i) = X_list(:, end);
end

% Calculate global truncation error
euler_truncation_error = abs(x_actual - euler_x_list);
midpoint_truncation_error = abs(x_actual - midpoint_x_list);

figure();
loglog(h_list, euler_truncation_error)
hold on;
loglog(h_list, midpoint_truncation_error)
xlabel('Step Size (h)');
ylabel('Global Truncation Error');
title('Global Truncation Error vs Step Size');
legend('Euler Error', 'Midpoint Error');
grid on;

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