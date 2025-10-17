clear;
clc;
clf;
close all;

%% Step Size Comparison
%%% forward_euler_fixed_step_integration
figure();
step_size_comparison(@forward_euler_fixed_step_integration)
title('Forward Euler H Step Comparison');
legend("Euler's Method H = 0.5", "Euler's Method H = 0.4", "Euler's Method H = 0.25",'Closed-Form Solution')

%%% backward_euler_fixed_step_integration
figure();
step_size_comparison(@backward_euler_fixed_step_integration)
title('Backward Euler H Step Comparison');
legend("Euler's Method H = 0.5", "Euler's Method H = 0.4", "Euler's Method H = 0.25",'Closed-Form Solution')

%%% explicit_midpoint_fixed_step_integration
figure();
step_size_comparison(@explicit_midpoint_fixed_step_integration)
title('Explicit Midpoint H Step Comparison');
legend("Explicit Midpoint H = 0.5", "Explicit Midpoint H = 0.4", "Explicit Midpoint H = 0.25",'Closed-Form Solution')

%%% implicit_midpoint_fixed_step_integration
figure();
step_size_comparison(@implicit_midpoint_fixed_step_integration);
title('Implicit Midpoint H Step Comparison');
legend("Implicit Midpoint H = 0.5", "Implicit Midpoint H = 0.4", "Implicit Midpoint H = 0.25",'Closed-Form Solution')

%% Explicit vs. Implicit Method Stability
%%% stability .38
h_ref = 0.38;
figure();
stability_comparison(h_ref);
title('.38 H Step Comparison');
legend("Forward Euler's Method H = 0.38","Backwards Euler's Method H = 0.38","Implicit Midpoint H = 0.38","Explicit Midpoint H = 0.38",'Closed-Form Solution')

%%% stability .45
h_ref = 0.45;
figure();
stability_comparison(h_ref);
title('.45 H Step Comparison');
legend("Forward Euler's Method H = 0.45","Backwards Euler's Method H = 0.45","Implicit Midpoint H = 0.45","Explicit Midpoint H = 0.45",'Closed-Form Solution')

%% Local Truncation Error Table
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, -1, num_points);

params = struct();
params.num_points = num_points;
params.t_ref = t_ref;
params.h_list = h_list;

step_funcs = {@forward_euler_step, @backward_euler_step, @explicit_midpoint_step, @implicit_midpoint_step};

%%% test_func01
test_1_local_slopes = zeros(1, length(step_funcs)+1);

% Calculate analytical difference
analytical_difference = zeros(1, length(h_list));
for i = 1:length(h_list)
    analytical_difference(i) = abs(solution01(t_ref + h_list(i)) - solution01(t_ref));
end
slope = polyfit(log(h_list), log(analytical_difference), 1);
test_1_local_slopes(1) = slope(1);

% Calculate local errors
for i = 1:length(step_funcs)
    error = local_truncation_error(step_funcs{i}, @rate_func01, @solution01, params);
    slope = polyfit(log(h_list), log(error), 1);
    test_1_local_slopes(i+1) = slope(1);
end

%%% test_func02
test_2_local_slopes = zeros(1, length(step_funcs)+1);

% Calculate analytical difference
analytical_difference = zeros(2, length(h_list));
for i = 1:length(h_list)
    analytical_difference(:, i) = abs(solution02(t_ref + h_list(i)) - solution02(t_ref));
end
analytical_difference = analytical_difference(1, :) + analytical_difference(2, :);
slope = polyfit(log(h_list), log(analytical_difference), 1);
test_2_local_slopes(1) = slope(1);

% Calculate local errors
for i = 1:length(step_funcs)
    error = local_truncation_error(step_funcs{i}, @rate_func02, @solution02, params);
    error = error(1, :) + error(2, :);
    slope = polyfit(log(h_list), log(error), 1);
    test_2_local_slopes(i+1) = slope(1);
end

%% Local Truncation Error Plot
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, 1, num_points);

params = struct();
params.num_points = num_points;
params.t_ref = t_ref;
params.h_list = h_list;

%%% Local Truncation Error Plot 1
euler_truncation_error = local_truncation_error(@forward_euler_step, @rate_func01, @solution01, params);
midpoint_truncation_error = local_truncation_error(@explicit_midpoint_step, @rate_func01, @solution01, params);

% Calculate analytical difference
analytical_difference = zeros(1, length(h_list));
for i = 1:length(h_list)
    analytical_difference(i) = abs(solution01(t_ref + h_list(i)) - solution01(t_ref));
end

% Create Fit Lines
p_analytical = polyfit(log(h_list), log(analytical_difference), 1);
y_hat_analytical = exp(p_analytical(1) * log(h_list) + p_analytical(2));

p_euler = polyfit(log(h_list), log(euler_truncation_error), 1);
y_hat_euler = exp(p_euler(1) * log(h_list) + p_euler(2));

p_midpoint = polyfit(log(h_list), log(midpoint_truncation_error), 1);
y_hat_midpoint = exp(p_midpoint(1) * log(h_list) + p_midpoint(2));

% Plot
figure();
loglog(h_list, y_hat_analytical, '-r')
hold on;
loglog(h_list, y_hat_euler, '-g')
loglog(h_list, y_hat_midpoint, '-b')

loglog(h_list, analytical_difference, '+r');
loglog(h_list, euler_truncation_error, '+g');
loglog(h_list, midpoint_truncation_error, '+b');
xlabel('Step Size (h)'); ylabel('Local Truncation Error'); title('Local Truncation Error for Explicit Methods'); 
legend('Analytical Difference', 'Euler Error', 'Midpoint Error');
grid on;

%%% Local Truncation Error Plot 2
backward_euler_truncation_error = local_truncation_error(@backward_euler_step, @rate_func01, @solution01, params);
implicit_midpoint_truncation_error = local_truncation_error(@implicit_midpoint_step, @rate_func01, @solution01, params);

figure();
loglog(h_list, euler_truncation_error)
hold on;
loglog(h_list, midpoint_truncation_error)
loglog(h_list, backward_euler_truncation_error);
loglog(h_list, implicit_midpoint_truncation_error);
xlabel('Step Size (h)'); ylabel('Local Truncation Error'); title('Local Truncation Error for All Methods'); legend('Forward Euler Error', 'Explicit Midpoint Error', 'Backward Euler Error', 'Implicit Midpoint Error');
grid on;

%% Global Truncation Error Table
num_points = 100;
h_list = logspace(-4, -1, num_points);
t_span = [0, 10];

params = struct();
params.num_points = num_points;
params.h_list = h_list;

step_funcs = {@forward_euler_step, @backward_euler_step, ...
              @explicit_midpoint_step, @implicit_midpoint_step};

num_methods = length(step_funcs);

test_1_global_slopes = zeros(1, num_methods);
test_2_global_slopes = zeros(1, num_methods);

%% Test Function 1 
for i = 1:num_methods
    [error, ~] = global_truncation_error(step_funcs{i}, @rate_func01, t_span, @solution01, params);
    slope = polyfit(log(h_list), log(error), 1);
    test_1_global_slopes(i) = slope(1);
end
%% Test Function 2 
for i = 1:num_methods
    [error, ~] = global_truncation_error(step_funcs{i}, @rate_func02, t_span, @solution02, params);
    slope = polyfit(log(h_list), log(error), 1);
    test_2_global_slopes(i) = slope(1);
end
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

%% Step Size Comparison Function
function step_size_comparison(fsi_func)
    t_span = [0, 6];
    h_list = [0.5, 0.4, 0.25];
    colors = ['b', 'k', 'm', 'r.'];
    X0 = solution01(0);

    for i = 1:length(h_list)
        [t_list, X_list, ~, ~] = fsi_func(@rate_func01, t_span, X0, h_list(i));
        plot(t_list, X_list, colors(i));
        hold on;
    end
    
    solution_X_list = zeros(length(X0), length(t_list));
    for i = 1:length(t_list)
        solution_X_list(i) = solution01(t_list(i));
    end
    plot(t_list, solution_X_list, colors(4));
    
    ylim([-2,2]);
    xlabel('Time');
    ylabel('Solution');
    grid on;
end

%% Stability Function
function stability_comparison(h_ref)
    fsi_funcs = {@forward_euler_fixed_step_integration, @backward_euler_fixed_step_integration, @implicit_midpoint_fixed_step_integration, @explicit_midpoint_fixed_step_integration};
    t_span = [0, 20];
    colors = ['b','g','r','m','k'];
    X0 = solution01(0);
    
    for i = 1:length(fsi_funcs)
        [t_list, X_list] = fsi_funcs{i}(@rate_func01, t_span, X0, h_ref);
        plot(t_list, X_list, colors(i));
        hold on;
    end

    solution_X_list = zeros(length(X0), length(t_list));
    for i = 1:length(t_list)
        solution_X_list(i) = solution01(t_list(i));
    end
    plot(t_list, solution_X_list, colors(5));

    ylim([-2,2])
    xlabel('Time');
    ylabel('Solution');
    grid on;
end

%% Local Truncation Error Function
function truncation_error = local_truncation_error(func, test_func, solution_func, params)
    %unpack values from struct (if fields in struct have been set)
    num_points = 100;
    if isfield(params,'num_points')
        num_points = params.num_points;
    end
    t_ref = 0.492;
    if isfield(params,'t_ref')
        t_ref = params.t_ref;
    end
    h_list = logspace(-5, 1, num_points);
    if isfield(params,'h_list')
        h_list = params.h_list;
    end
    
    % Calculate
    X0 = solution_func(t_ref);
    x_list = zeros(length(X0), length(h_list));
    x_actual_list = zeros(length(X0), length(h_list));
    
    for i = 1:num_points
        x_list(:, i) = func(test_func, t_ref, X0, h_list(i));
        x_actual_list(:, i) = solution_func(t_ref + h_list(i));
    end
    truncation_error = abs(x_actual_list - x_list);
end

%% Global Truncation Error Function (for rate_func01)
function [truncation_error, h_avg_list] = global_truncation_error(step_func, test_func, t_span, solution_func, params)

    num_points = 100;
    if isfield(params, 'num_points')
        num_points = params.num_points;
    end

    h_list = logspace(-5, 1, num_points);
    if isfield(params, 'h_list')
        h_list = params.h_list;
    end

    X0 = solution_func(t_span(1));
    dim = length(X0); 

    x_approx_list = zeros(dim, length(h_list));
    x_actual_list = zeros(dim, length(h_list));
    h_avg_list = zeros(1, length(h_list));

    for i = 1:length(h_list)
        [t_list, X_list, h_avg, ~] = fixed_step_integration(test_func, step_func, t_span, X0, h_list(i));
        x_approx_list(:, i) = X_list(:, end);
        x_actual_list(:, i) = solution_func(t_span(2));
        h_avg_list(i) = h_avg;
    end

    truncation_error = vecnorm(x_actual_list - x_approx_list, 2, 1);
end