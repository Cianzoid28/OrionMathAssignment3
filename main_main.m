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
h_list = logspace(-5, 1, num_points);

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
slope = polyfit([log(h_list); log(h_list)], log(analytical_difference), 1);
test_2_local_slopes(1) = slope(1);

% Calculate local errors
for i = 1:length(step_funcs)
    error = local_truncation_error(step_funcs{i}, @rate_func02, @solution02, params);
    slope = polyfit([log(h_list), log(h_list)], log(error), 1);
    test_2_local_slopes(i+1) = slope(1);
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