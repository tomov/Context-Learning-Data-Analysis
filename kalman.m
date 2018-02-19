function [x, P] = kalman(z, x, P, F, B, u, Q, R)

%
% x_t = F * x_t-1 + B u + w_t
% z_t = H * x_t + v_t
%
% w_t ~ N(0, Q)
% v_t ~ N(0, R)
%
% x_t = state at time t
% z_t = observation/measurement at time t
% F = state transition model (maps x to next x)
% B = control input model (maps u to x)
% u = control input
% w_t = process noise
% R = process noise covariance
% v_t = observation/measurement noise
% Q = measurement noise covariance
% H = observation/measurement model (maps x to z)
%

% Prediction
% Note x and P now correspond to x_hat_t|t-1 and P_t|t-1
%
x = F * x + B * u;
P = F * P * F' + Q;

% Observation update
% Note x and P now correspond to x_hat_t|t and P_t|t
%
K = P * H' * inv(H * P * H' + R);
x = x + K * (z - H * x);
P = P - K * H * P;


