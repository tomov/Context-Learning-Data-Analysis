function [x, P, z_mean, z_var] = kalman(z, x, P, F, B, u, Q, H, R)

% Generic Kalman filter update from time t-1 to time t.
% Following notation from Wikipedia.
%
% INPUT:
% z = z_t = observation at time t
% x = x_hat_t-1|t-1 = (posterior) mean of state x at time t-1 = E[x_t-1 | z_1..t-1]
% P = P_t-1|t-1 = (posterior) covariance of state x at time t-1 = Var[x_t-1 | z_1..t-1]
% F = F_t = state transition model at time t
% B = B_t = control input model at time t
% u = u_t = control input at time t
% Q = Q_t = covariance of process noise at time t
% H = H_t = observation model at time t
% R = R_t = covariance of observation noise at time t
%
% OUTPUT:
% x = x_hat_t|t = (posterior) mean of state x at time t = E[x_t | z_1..t]
% P = P_t|t = (posterior) covariance of state x at time t = Var[x_t | z_1..t]
% z_mean = (posterior) expected observation at time t = E[z_t | x_hat_t|t]
% z_var = (posterior) covariance of observation at time t = Var[z_t | x_hat_t|t]

%
% x_t = F_t * x_t-1 + B_t u_t + w_t
% z_t = H_t * x_t + v_t
%
% w_t ~ N(0, Q_t)
% v_t ~ N(0, R_t)
%
% x_t = state at time t
% z_t = observation/measurement at time t
% F_t = state transition model (maps x to next x)
% B_t = control input model (maps u to x)
% u_t = control input
% w_t = process noise
% R_t = process noise covariance
% v_t = observation/measurement noise
% Q_t = measurement noise covariance
% H_t = observation/measurement model (maps x to z)
%

% Prediction
% Note x and P now correspond to x_hat_t|t-1 and P_t|t-1
%
z_mean = H * x;

x = F * x + B * u;
P = F * P * F' + Q;

z_var = H * P * H' + R; % TODO SAM why use this variance?

% Observation update
% Note x and P now correspond to x_hat_t|t and P_t|t
%
K = P * H' / (H * P * H' + R);
x = x + K * (z - H * x);
P = P - K * H * P;

