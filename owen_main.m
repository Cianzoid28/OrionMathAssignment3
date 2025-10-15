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

%% backward_euler_fixed_step_integration
t_span = [0, 6];
h_list = [0.5, 0.4, 0.25];
colors = ['b', 'k', 'm', 'r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = backward_euler_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
title('Backward Euler H Step Comparison');
legend("Euler's Method H = 0.5", "Euler's Method H = 0.4", "Euler's Method H = 0.25",'Closed-Form Solution')
grid on;

%% explicit_midpoint_fixed_step_integration
t_span = [0, 6];
h_list = [0.5, 0.4, 0.25];
colors = ['b', 'k', 'm', 'r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
title('Explicit Midpoint H Step Comparison');
legend("Explicit Midpoint H = 0.5", "Explicit Midpoint H = 0.4", "Explicit Midpoint H = 0.25",'Closed-Form Solution')
grid on;

%% Implicit_midpoint_fixed_step_integration
t_span = [0, 6];
h_list = [0.5, 0.4, 0.25];
colors = ['b', 'k', 'm', 'r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = implicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
title('Implicit Midpoint H Step Comparison');
legend("Implicit Midpoint H = 0.5", "Implicit Midpoint H = 0.4", "Implicit Midpoint H = 0.25",'Closed-Form Solution')
grid on;
%% forward_euler_fixed_step_stability
t_span = [0, 20];
h_list = [0.38];
colors = ['b','r.'];
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
plot(t_list, solution_X_list, colors(2))
title('Forward Euler H Step Comparison');
legend("Euler's Method H = 0.38",'Closed-Form Solution')
grid on;

%% backward_euler_fixed_step_stability
t_span = [0, 20];
h_list = [0.38];
colors = ['b','r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = backward_euler_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
plot(t_list, solution_X_list, colors(2))
title('Backward Euler H Step Comparison');
legend("Euler's Method H = 0.38",'Closed-Form Solution')
grid on;

%% explicit_midpoint_fixed_step_stability
t_span = [0, 20];
h_list = [0.38];
colors = ['b', 'r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
plot(t_list, solution_X_list, colors(2))
title('Explicit Midpoint H Step Comparison');
legend("Explicit Midpoint H = 0.38",'Closed-Form Solution')
grid on;

%% Implicit_midpoint_fixed_step_stability
t_span = [0, 20];
h_list = [0.38];
colors = ['b','r.'];
% rate_func01
X0 = solution01(0);
figure();
for j = 1:length(h_list)

    [t_list, X_list, h_avg, num_evals] = implicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_list(j));
    
    
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
plot(t_list, solution_X_list, colors(2))
title('Implicit Midpoint H Step Comparison');
legend("Implicit Midpoint H = 0.38",'Closed-Form Solution')
grid on;

%% %% Defined functions
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