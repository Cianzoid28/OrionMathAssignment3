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


%% p
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, -1, num_points);

euler_truncation_error = local_truncation_error(@forward_euler_step, num_points, t_ref, h_list);
midpoint_truncation_error = local_truncation_error(@explicit_midpoint_step, num_points, t_ref, h_list);
backward_euler_truncation_error = local_truncation_error(@backward_euler_step, num_points, t_ref, h_list);
implicit_midpoint_truncation_error = local_truncation_error(@implicit_midpoint_step, num_points, t_ref, h_list);

% Compute logs
log_h_list = log(h_list);
log_euler_error_list = log(euler_truncation_error);
log_midpoint_error_list = log(midpoint_truncation_error);
log_backward_error_list = log(backward_euler_truncation_error);
log_implicit_midpoint_error_list = log(implicit_midpoint_truncation_error);


% Linear fit
coeffsFE = polyfit(log_h_list, log_euler_error_list, 1);
coeffsFME = polyfit(log_h_list, log_midpoint_error_list, 1);
coeffsBE = polyfit(log_h_list, log_backward_error_list, 1);
coeffsBME = polyfit(log_h_list, log_implicit_midpoint_error_list, 1);

% Extract slope (p)
pFE = coeffsFE(1);
pMFE = coeffsFME(1);
pBE = coeffsBE(1);
pBME = coeffsBME(1);



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

%% Fixed Step Integration
function [X_list, solution_X_list, t_list] = calculate_fixed_step_integration(step_func, test_func, solution_func, t_span, h)
    X0 = solution_func(0);
    [t_list, X_list, h_avg, num_evals] = step_func(test_func, t_span, X0, h);
    solution_X_list = zeros(length(X0), length(t_list));
    for i = 1:length(t_list)
        solution_X_list(:, i) = solution_func(t_list(i));
    end
end

%% Local Truncation Error Function (for rate_func01)
function truncation_error = local_truncation_error(func, num_points, t_ref, h_list)
    X0 = solution01(t_ref);
    x_list = zeros(1, length(h_list));
    x_actual_list = zeros(1, length(h_list));
    for i = 1:num_points
        [X1, num_evals] = func(@rate_func01, t_ref, X0, h_list(i));
        x_list(i) = X1;
        x_actual_list(i) = solution01(t_ref + h_list(i));
    end
    truncation_error = abs(x_actual_list - x_list);
end

%% Global Truncation Error Function (for rate_func01)
function [truncation_error, h_avg_list] = global_truncation_error(func, func2, num_points, t_span, h_list)
    X0 = solution01(t_span(1));
    x_actual = solution01(t_span(length(t_span)));
    x_list = zeros(1, length(h_list));
    h_avg_list = zeros(1, length(h_list));
    if isa(func2, 'function_handle')
        for i = 1:num_points
            [t_list, X_list, h_avg, num_evals] = func(@rate_func01, func2, t_span, X0, h_list(i));
            x_list(i) = X_list(:, end);
            h_avg_list(i) = h_avg;
        end
    else
        for i = 1:num_points
            [t_list, X_list, h_avg, num_evals] = func(@rate_func01, t_span, X0, h_list(i));
            x_list(i) = X_list(:, end);
            h_avg_list(i) = h_avg;
        end
    end
    truncation_error = abs(x_actual - x_list);
end

%% Local Truncation Error Function (for rate_func02)
function truncation_error = local_truncation_error2(func, num_points, t_ref, h_list)
    X0 = solution02(t_ref);
    x_list = zeros(1, length(h_list));
    x_actual_list = zeros(1, length(h_list));
    for i = 1:num_points
        [X1, num_evals] = func(@rate_func01, t_ref, X0, h_list(i));
        x_list(i) = X1;
        x_actual_list(i) = solution02(t_ref + h_list(i));
    end
    truncation_error = abs(x_actual_list - x_list);
end
