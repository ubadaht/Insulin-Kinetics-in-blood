clear all; clc; close all;
global kin k12 k21 k13 k31

kin = 2;
k12 = 1;
k21 = 1.5;
k13 = 2;
k31 = 1.75;

x0 = [0,0,0,25]; %initial conditions
tspan = [0,10]; %time span upto 10 hrs

% system of diff equations
Fx = @(t,x) [kin*x(4) + k31*x(3) + k21*x(2) - k13*x(1) - k12*x(1);
             k12*x(1) - k21*x(2);
             k13*x(1) - k31*x(3);
             -kin*x(4)];

h = 0.1;
[t_exp,X_exp] = euler_exp(Fx,tspan,x0,h);
% Define tolerance
tolerance = 1e-6;
% specify tolerance
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
[t,X] = ode45(Fx,tspan,x0,options);
[t_heun,X_heun] = heun(Fx,tspan,x0,h);
[t_imp,X_imp]=eul_imp(Fx,@jacf,tspan,x0,h);

%plots
figure()
plot(t_exp,X_exp(:,1),t_heun,X_heun(:,1),t_imp,X_imp(:,1),t,X(:,1))
legend('eul-exp','heun','eul-imp','ode45')
xlabel('time(hours)')
ylabel('Insulin (U units)')
title('Blood')
grid on;


figure()
plot(t_exp,X_exp(:,2),t_heun,X_heun(:,2),t_imp,X_imp(:,2),t,X(:,2))
legend('eul-exp','heun','eul-imp','ode45')
xlabel('time(hours)')
ylabel('Insulin (U units)')
title('Kidneys')
grid on;

figure()
plot(t_exp,X_exp(:,3),t_heun,X_heun(:,3),t_imp,X_imp(:,3),t,X(:,3))
legend('eul-exp','heun','eul-imp','ode45')
xlabel('time(hours)')
ylabel('Insulin (U units)')
title('Pancreas')
grid on;

figure()
plot(t_exp,X_exp(:,4),t_heun,X_heun(:,4),t_imp,X_imp(:,4),t,X(:,4))
legend('eul-exp','heun','eul-imp','ode45')
xlabel('time(hours)')
ylabel('Insulin (U units)')
title('Abdomen')
grid on;

% Extract last values for each compartment
ode45_last_values = X(end,:);
exp_last_values = X_exp(end,:);
heun_last_values = X_heun(end,:);
imp_last_values = X_imp(end,:);

% Steady state value for blood compartment
steady_state_blood = [ode45_last_values(1), exp_last_values(1), heun_last_values(1), imp_last_values(1)];
disp(['Steady state value for blood compartment: ', num2str(steady_state_blood)]);


% Calculate relative errors
relative_error_exp = abs((exp_last_values - ode45_last_values) ./ ode45_last_values)*100;
relative_error_heun = abs((heun_last_values - ode45_last_values) ./ ode45_last_values)*100;
relative_error_imp = abs((imp_last_values - ode45_last_values) ./ ode45_last_values)*100;

% Display relative errors
disp('Relative error at last index for each compartment:');
disp('Compartment:   Euler Explicit   Heun    Euler Implicit');
disp(['Blood:         ', num2str(relative_error_exp(1)), '    ', num2str(relative_error_heun(1)), '    ', num2str(relative_error_imp(1))]);
disp(['Kidneys:       ', num2str(relative_error_exp(2)), '    ', num2str(relative_error_heun(2)), '    ', num2str(relative_error_imp(2))]);
disp(['Pancreas:      ', num2str(relative_error_exp(3)), '    ', num2str(relative_error_heun(3)), '    ', num2str(relative_error_imp(3))]);
disp(['Abdomen:       ', num2str(relative_error_exp(4)), '         ', num2str(relative_error_heun(4)), '     ', num2str(relative_error_imp(4))]);


% Arrange data for bar plot
relative_errors = [relative_error_exp; relative_error_heun; relative_error_imp];
method_names = {'Euler Exp', 'Heun', 'Euler Imp'};
compartment_names = {'Blood', 'Kidneys', 'Pancreas', 'Abdomen'};

% Create separate bar plots for each compartment using subplots
figure();
for i = 1:4
    subplot(2, 2, i);
    bar(relative_errors(:, i));
    title(compartment_names{i});
    xlabel('Method');
    ylabel('Relative Error');
    set(gca, 'XTickLabel', method_names);
end
sgtitle('Relative Error at Last Index for Each Compartment');