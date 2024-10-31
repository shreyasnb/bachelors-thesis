n = 100;
h = 1e-1;
%% Initialize states
I = eye(3);
O = zeros(3,3);

% FL matrices
A = [[O I]; [O O]];
B = [O; I];
kp = 10;
ki = 5;
K = [I*ki I*kp];

% Angular velocities
p=0; q=0; r=0;
Omega_0 = [p; q; r];

% Initial rotation matrix
r11 = cos(pi/2); r12 = 0; r13 = -sin(pi/2); 
r21 = 0; r22=1; r23=0; 
r31 = sin(pi/2); r32=0; r33=cos(pi/2);

%% Define states in [SO(3) x R^3]
R_0 = [[r11 r12 r13]; [r21 r22 r23]; [r31 r32 r33]];
x = zeros(n, 12, 1);
x_k = zeros(n, 12, 1);
x(1, :, :) = [columnize(R_0); Omega_0];
x_k(1, :, :) = [columnize(R_0); Omega_0];
z = zeros(n, 6, 1);
z(1,:,:) = [alpha(R_0); Omega_0];
g = [zeros(9,3); I];

%% Nonlinear dynamics on ODE45
tspan = [0 n*h];
[t,y] = ode45(@(t, x) odefun(t, x), tspan, x(1,:,:));

R = zeros(length(y), 3, 3);
err = zeros(length(y),1);
for i=1:length(y)
    R(i, :, :) = matrix(y(i, 1:9)');
    err(i,1) = trace(I - matrix(y(i, 1:9)'));
end

%% FL dynamics using forward Euler scheme
for k=1:n-1
    z(k+1,:,:) = ([[I O]; [O I]] +h*(A-B*K))*(z(k,:,:)');
    x_k(k+1,:,:) = [columnize(expm(skew(z(k+1, 1:3,1))))', z(k+1,4:6,:)];
end

R_k = zeros(n,3,3);
err_k = zeros(n,1);
for i=1:n
    R_k(i,:,:) = matrix(x_k(i,1:9,1));
    err_k(i,1) = trace(I-matrix(x_k(i,1:9,1)));
end


figure
plot(t, err, 'DisplayName','ODE45');
hold on
grid on
plot(linspace(0,n*h,n), err_k, 'DisplayName', 'Forward Euler on FL system')
xlabel('Time (s)', 'Interpreter','latex')
ylabel('tr(I-R)', 'Interpreter','latex')
legend('Interpreter','latex')

figure
hold on
grid on
plot(t, y(:,10), 'DisplayName','$\Omega_1$')
plot(t, y(:,11), 'DisplayName','$\Omega_2$')
plot(t, y(:,12), 'DisplayName','$\Omega_3$')
xlabel('Time (s)', 'Interpreter','latex')
ylabel('Angular velocities', 'Interpreter','latex')
legend('Interpreter','latex')

function dxdt = odefun(t, x)
    %% Choose PI Controller Gains
    kp = 10;
    ki = 5;
    g = [zeros(9,3); eye(3)];
    dxdt = f(x) + g*(-ki*alpha(x)-kp*x(10:12));
end

function uncap = vec(M)
    uncap = [M(3,2); M(1,3); M(2,1)];
end

function mat = matrix(v)
    mat = zeros(3,3);
    k=1;
    for i=1:3
        for j=1:3
            mat(j,i) = v(k);
            k=k+1;
        end
    end
end

function v = alpha(x)
    R = matrix(x(1:9));
    v = vec(logm(R));
end

function cap = skew(w)
    cap = [[0 -w(3) w(2)]; [w(3) 0 -w(1)]; [-w(2) w(1) 0]];
end

function sys = f(x)
    p = x(10);
    q = x(11);
    r = x(12);
    sys = zeros(12,1);
    sys(1:3,1) = r*x(4:6) - q*x(7:9);
    sys(4:6,1) = p*x(7:9) - r*x(1:3);
    sys(7:9,1) = q*x(1:3) - p*x(4:6);
end

function vec = columnize(M)
    vec = zeros(9, 1);
    k=1;
    for i=1:3
        for j=1:3
            vec(k,1) = M(j,i);
            k=k+1;
        end
    end
end