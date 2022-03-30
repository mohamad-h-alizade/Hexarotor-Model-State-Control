clear; close all; clc;
syms phi theta psi z z_dot P Q R u1 u2 u3 u4 Ixx Iyy Izz Ixy Ixz Iyz Ip b L d m g

I = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
x_ss = zeros(8, 1);
u_ss = [m*g/b; 0; 0; 0];
pm_0 = [0.01; -0.01; 0.01];
eul_0 = [0.1; -0.1; 0.1];
rdot_0 = [0; 0; 0.01];
r_0 = [0; 0; 0.1];

w = [P;Q;R];
u = [u1; u2; u3; u4];

x = [z; z_dot; phi; theta; psi; P; Q; R];
multi = [0.1667 0 0.3333 -0.1667; 0.1667 -0.2887 0.1667 0.1667; 0.1667 -0.2887 -0.1667 -0.1667;...
    0.1667 0 -0.3333 0.1667; 0.1667 0.2887 -0.1667 -0.1667; 0.1667 0.2887 0.1667 0.1667];
W = sqrt(multi*u);
gu = W(1) - W(2) + W(3) - W(4) + W(5) - W(6);
h = [0;0;Ip*gu];

x_dot(1) = x(2);
x_dot(2) = g - cos(phi)*cos(theta)*(b/m)*u1;
X = 1 / cos(theta) * [cos(theta) sin(theta)*sin(phi) sin(theta)*cos(phi); 0 cos(theta)*cos(phi) -cos(theta)*sin(phi); 0 sin(phi) cos(phi)] * [P;Q;R];
x_dot(3) = X(1);
x_dot(4) = X(2);
x_dot(5) = X(3);
X = pinv(I) * ([L*b*u2; L*b*u3; d*u4] - cross(w, (I*w+h)));
x_dot(6) = X(1);
x_dot(7) = X(2);
x_dot(8) = X(3);

%% linearization
A_sym = jacobian(x_dot, x);
B_sym = jacobian(x_dot, u);
params = [x; u; Ixx; Iyy; Izz; Ixy; Ixz; Iyz; Ip; b; L; d; m; g];
vals = [x_ss; u_ss; 4.85e-3; 4.95e-3; 8.81e-3; 0; 0; 0; 3.36e-5; 2.92e-6; 0.2; 1.12e-7; 0.5; 9.81];
A = double(subs(A_sym, params, vals));
B = double(subs(B_sym, params, vals));

C = [1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0];
D = zeros(4, 4);
%% Controller
p_1 = [-1;
    -5;
    -1;
    -1;
    -1;
    -5;
    -5;
    -5];
k1 = place(A,B,p_1)
req = [5;
    0;
    0;
    -1];
r = [10;
    pi/6;
    -pi/6;
    10*pi/18];
k2 = inv(C*inv(-A+B*k1)*B)
p_2 = [-5;-10;-5;-5;-5;-10;-10;-10];
L_observer = (place(A',C',p_2))'

%% Integrator coefficients
A_bar = [A zeros(8,4);-C zeros(4,4)];
B_bar = [B;zeros(4,4)];
p_3 = [-0.9;-0.8;-0.9;-0.9;-0.9;-0.8;-0.8;-0.8;-1;-1;-1;-1];
k_bar = place(A_bar,B_bar,p_3);
k3 = k_bar(:,1:8)
k4 = k_bar(:,9:12)

%% parameters for simulink model
g = 9.81;
m = 0.5;
L = 0.2;
b = 2.92e-6;
d = 1.12e-7;
I_xx = 4.85e-3;
I_yy = 4.95e-3;
I_zz = 8.81e-3;
I_xy = 0;
I_xz = 0;
I_yz = 0;
I_p = 3.36e-5;
I = [I_xx I_xy I_xz;
    I_xy I_yy I_yz;
    I_xz I_yz I_zz];
%% linear model simulation
out = sim("Linear.slx");

figure;
subplot(2, 2, 1)
plot(out.tout, out.z)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi)
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Linear response')
%% non linear model simulation (tau = 0.1)
tau = 0.1;
out = sim("Hex.slx");

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 0.1s')

%% non linear model simulation (tau = 1)
tau = 1;
out = sim("Hex.slx");

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 1s')

%% non linear model simulation (sinusoidal tracking with tau = 0.1)
tau = 0.1;

A_z = 5; 
A_phi = 0; 
A_theta = 0; 
A_psi = 0;
out = sim("Sinus.slx", 600);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('time(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 0.1s and sinusoidal input')
%%
A_z = 0; 
A_phi = 1; 
A_theta = 0; 
A_psi = 0;
out = sim("Sinus.slx", 600);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 0.1s and sinusoidal input')
%%
A_z = 0; 
A_phi = 0; 
A_theta = 1.3; 
A_psi = 0;
out = sim("Sinus.slx", 600);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 0.1s and sinusoidal input')
%%
tau=0.1;

A_z = 0; 
A_phi = 0; 
A_theta = 0; 
A_psi = 0.07;
out = sim("Sinus.slx", 600);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 0.1s and sinusoidal input')

%%
tau = 1;

A_z = 0;
A_phi = 0.2;
A_theta = 0;
A_psi = 0;
out = sim("Sinus.slx", 30);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('time(s)')
ylabel('z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('time(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 1s and sinusoidal input')

A_z = 0; 
A_phi = 0;
A_theta = 0.2;
A_psi = 0;
out = sim("Sinus.slx", 30);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 1s and sinusoidal input')

A_z = 0; 
A_phi = 0; 
A_theta = 0; 
A_psi = 1;
out = sim("Sinus.slx", 30);

figure;
subplot(2, 2, 1)
plot(out.tout, out.z1)
hold on
plot(out.tout, out.zd)
xlabel('t(s)')
ylabel('Z(m)')

subplot(2, 2, 2)
plot(out.tout, out.phi1)
hold on
plot(out.tout, out.phid)
legend('Actual', 'Desired')
xlabel('t(s)')
ylabel('\phi(deg)')

subplot(2, 2, 3)
plot(out.tout, out.theta1)
hold on
plot(out.tout, out.thetad)
xlabel('t(s)')
ylabel('\theta(deg)')

subplot(2, 2, 4)
plot(out.tout, out.psi1)
hold on
plot(out.tout, out.psid)
xlabel('t(s)')
ylabel('\psi(deg)')

sgtitle('Non-linear response with \tau = 1s and sinusoidal input')
