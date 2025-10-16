clear;
clc;
clf;
close all;

%% forward_euler_step
t = 0;
h = 0.3;

% rate_func01
X0 = solution01(t);
labels = ["Forward Euler Step rate_func01:", "Expected Step rate_func01:"];
results = [forward_euler_step(@rate_func01, t, X0, h), solution01(t+h)];
disp([labels; results]);

% rate_func02
X0 = solution02(t);
labels = ["Forward Euler Step rate_func02:", "Expected Step rate_func02:"];
results = [forward_euler_step(@rate_func02, t, X0, h), solution02(t+h)];
disp([labels; results]);

%% forward_euler_fixed_step_integration
t_span = [0, 6];

% rate_func01
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@forward_euler_fixed_step_integration, @rate_func01, @solution01, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list);
xlabel('Time'); ylabel('Solution'); title('Test01 Forward Euler Integration Results'); legend("Euler's Method", 'Closed-Form Solution')
grid on;

% rate_func02
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@forward_euler_fixed_step_integration, @rate_func02, @solution02, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list);
xlabel('Time'); ylabel('Solution'); title('Test02 Forward Euler Integration Results');
grid on;

%% explicit_midpoint_fixed_step_integration

% rate_func01
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@explicit_midpoint_fixed_step_integration, @rate_func01, @solution01, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list);
xlabel('Time'); ylabel('Solution'); title('Test01 Explicit Midpoint Integration Results'); legend("Midpoint Method", 'Closed-Form Solution')
grid on; 

% rate_func02
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@explicit_midpoint_fixed_step_integration, @rate_func02, @solution02, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list);
xlabel('Time'); ylabel('Solution'); title('Test02 Explicit Midpoint Integration Results');
grid on;

%% Local Truncation Error
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, 1, num_points);

euler_truncation_error = local_truncation_error(@forward_euler_step, num_points, t_ref, h_list);
midpoint_truncation_error = local_truncation_error(@explicit_midpoint_step, num_points, t_ref, h_list);

figure();
loglog(h_list, euler_truncation_error)
hold on;
loglog(h_list, midpoint_truncation_error)
xlabel('Step Size (h)'); ylabel('Local Truncation Error'); title('Local Truncation Error vs Step Size'); legend('Euler Error', 'Midpoint Error');
grid on;

%% Global Truncation Error
t_span = [0, t_ref];

% Calculate global truncation error
[euler_truncation_error, h_avg_list] = global_truncation_error(@forward_euler_fixed_step_integration, 0, num_points, t_span, h_list);
midpoint_truncation_error = global_truncation_error(@explicit_midpoint_fixed_step_integration, 0, num_points, t_span, h_list);

figure();
loglog(h_avg_list, euler_truncation_error)
hold on;
loglog(h_avg_list, midpoint_truncation_error)
xlabel('Step Size (h)'); ylabel('Global Truncation Error'); title('Global Truncation Error vs Step Size'); legend('Euler Error', 'Midpoint Error');
grid on;

%% backward_euler_step
t = 0;
h = 0.3;

% rate_func01
X0 = solution01(t);
labels = ["Backwards Euler Step rate_func01:", "Expected Step rate_func01:"];
results = [backward_euler_step(@rate_func01, t, X0, h), solution01(t+h)];
disp([labels; results]);

% rate_func02
X0 = solution02(t);
labels = ["Backwards Euler Step rate_func02:", "Expected Step rate_func02:"];
results = [backward_euler_step(@rate_func02, t, X0, h), solution02(t+h)];
disp([labels; results]);

%% fixed_step_integration
t_span = [0, 6];

% backward_euler_step rate_func01
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@backward_euler_fixed_step_integration, @rate_func01, @solution01, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time'); ylabel('Solution'); title('Test01 Backward Euler Integration Results'); legend("Backward Euler's Method", 'Closed-Form Solution')
grid on;

% backward_euler_step rate_func02
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@backward_euler_fixed_step_integration, @rate_func02, @solution02, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time'); ylabel('Solution'); title('Test02 Backward Euler Integration Results');
grid on;

% implicit_midpoint_step rate_func01
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@implicit_midpoint_fixed_step_integration, @rate_func01, @solution01, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list)
xlabel('Time'); ylabel('Solution'); title('Test01 Implicit Midpoint Integration Results'); legend("Implicit Midpoint Method", 'Closed-Form Solution')
grid on;

% implicit_midpoint_step rate_func02
[X_list, solution_X_list, t_list] = calculate_fixed_step_integration(@implicit_midpoint_fixed_step_integration, @rate_func02, @solution02, t_span, h);

figure();
hold on;
plot(t_list, X_list);
plot(t_list, solution_X_list);
xlabel('Time'); ylabel('Solution'); title('Test02 Implicit Midpoint Integration Results');
grid on;

%% Local Truncation Error
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, 1, num_points);

% Calculate local truncation error
backward_euler_truncation_error = local_truncation_error(@backward_euler_step, num_points, t_ref, h_list);
implicit_midpoint_truncation_error = local_truncation_error(@implicit_midpoint_step, num_points, t_ref, h_list);

figure();
loglog(h_list, backward_euler_truncation_error);
hold on;
loglog(h_list, implicit_midpoint_truncation_error);
xlabel('Step Size (h)'); ylabel('Local Truncation Error'); title('Local Truncation Error vs Step Size'); legend('Backward Euler Error', 'Implicit Midpoint Error');
grid on;

%% Global Truncation Error
t_span = [0, t_ref];

% Calculate global truncation error
[euler_truncation_error, h_avg_list] = global_truncation_error(@fixed_step_integration, @forward_euler_step, num_points, t_span, h_list);
midpoint_truncation_error = global_truncation_error(@fixed_step_integration, @explicit_midpoint_step, num_points, t_span, h_list);
backward_euler_truncation_error = global_truncation_error(@fixed_step_integration, @backward_euler_step, num_points, t_span, h_list);
implicit_euler_truncation_error = global_truncation_error(@fixed_step_integration, @implicit_midpoint_step, num_points, t_span, h_list);

figure();
loglog(h_avg_list, euler_truncation_error);
hold on;
loglog(h_avg_list, midpoint_truncation_error);
loglog(h_avg_list, backward_euler_truncation_error);
loglog(h_avg_list, implicit_midpoint_truncation_error);
xlabel('Step Size (h)'); ylabel('Global Truncation Error'); title('Global Truncation Error vs Step Size'); legend('Forward Euler Error', 'Explicit Midpoint Error', 'Backward Euler Error', 'Implicit Midpoint Error');
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