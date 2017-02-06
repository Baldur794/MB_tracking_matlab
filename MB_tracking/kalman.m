% omega_s = 1;
% zeta_s = 0.1;
% Noise intensity
pros_noise_acc = 0.5;
meas_noise_pos = 0.5;
v1 = [pros_noise_acc pros_noise_acc]';
v2 = [meas_noise_pos meas_noise_pos]';

% System
A = zeros(4);
A(1,2) = 1;
% A(2,5) = omega_s;
A(3,4) = 1;
% A(4,7) = omega_s;
% A(5,6) = 1;
% A(6,5) = -omega_s^2;
% A(6,6) = -2*zeta_s*omega_s;
% A(7,8) = 1;
% A(8,7) = -omega_s^2;
% A(8,8) = -2*zeta_s*omega_s;


B = zeros(4,1);
Bv = zeros(4,2);
Bv(2,1) = 1;
Bv(4,2) = 1;
C = zeros(2,4);
C(1,1) = 1;
C(2,3) = 1;
D = 0;

V1 = diag(v1);
V2 = diag(v2);

% System init
sys_init = [0 0 0 0];


% Parameters
STEP_TIME = 1e-4;
SIM_TIME = 100;
SEED_VEC_v1 = [1 2];
SEED_VEC_v2 = [3 4];
T_s = 1/50;
f_max = 50;
T_c = STEP_TIME;%1/50*6/f_max;

v1_switch = 1;
v2_switch = 1;

systemus_c = ss(A,[B Bv],C,D);
%[~,L,P] = kalman(systemus_c,V1,V2);

T_s = 1;
systemus_d = c2d(systemus_c,T_s);
[~,L,P] = kalman(systemus_d,V1,V2);
F = systemus_d.a;
G = systemus_d.b(:,size(B,2));
Gv = systemus_d.b(:,end-size(Bv,2)+1:end);

%% Predict states
x_old = zeros(4,1);
y_old = zeros(2,1);
x_new = F*x_old+F*L*(y_old-C*x_old)


%% Plotting
time = simout.Time;
pos = simout.Data(:,1:4);
vel = simout.Data(:,5:8);

figure();
plot(pos(:,1),pos(:,3));
hold on
plot(pos(:,2),pos(:,4));

%%
sim('kalman_sin');
