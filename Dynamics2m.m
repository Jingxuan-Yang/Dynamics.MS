%Author:  JingXuan Yang
%E-mail:  yangjingxuan@stu.hit.edu.cn
%Date:    2019.03.29
%Project: Dynamics of Mechanical System Homework 2
%Purpose: Dynamics of 2 DOF manipulator
%Note   : All angles in this script are in radius

clear;  
clc;

%initial data
%problem 1
M11 = 10;
M21 = 0;

%problem 2
M12 = 0;
M22 = 10;

%problem 3
%number given in this problem
% M13d = [5, 30, 30, -4, -6, -25.05, -17.64,  0];
% M23d = [2.29, 13.88, 14.38, 0.4926, 2.628,  -10.28, -9.423, 0];

%construct new torque from debugging
M13d = [5, 30, 29.95, -4.533, -6.8, -22.455, -15.04, 0];
M23d = [2.28, 13.96, 14.114, 0.4986, 2.285, -8.8097, -7.323, 0];
theta1 = [0, pi/80, 4*pi/80, pi/8,3*pi/8, 36*pi/80, 39*pi/80, pi/2];

%number of loop time
%for problem 1 and 2
n = 333;
%for problem 3
m = 101;

%step for Runge-Kutta Method
h = 0.03;

%initialize the matrices
theta11 = zeros(1,n);
theta21 = zeros(1,n);
dtheta11 = zeros(1,n);
dtheta21 = zeros(1,n);
ddtheta11 = zeros(1,n);
ddtheta21 = zeros(1,n);

theta12 = zeros(1,n);
theta22 = zeros(1,n);
dtheta12 = zeros(1,n);
dtheta22 = zeros(1,n);
ddtheta12 = zeros(1,n);
ddtheta22 = zeros(1,n);

theta13 = zeros(1,m);
theta23 = zeros(1,m);
dtheta13 = zeros(1,m);
dtheta23 = zeros(1,m);
ddtheta13 = zeros(1,m);
ddtheta23 = zeros(1,m);

Theta1 = zeros(n,4);
Theta2 = zeros(n,4);
Theta3 = zeros(m,4);

%calculating
%time for problem 1 and 2
t = h*(0:1:n-1);

%solving problem 1 and 2
for i = 1:n
    
    if i == 1
        %initialize the matrices
        Theta1(1,:) = zeros(1,4);        
        Theta2(1,:) = zeros(1,4);
    else
        Theta1(i,:) = RK(Theta1(i-1,:),M11,M21,h);
        theta11(i) = Theta1(i,1);
        theta21(i) = Theta1(i,2);
        dtheta11(i) = Theta1(i,3);
        dtheta21(i) = Theta1(i,4);
        ddtheta11(i) = f1(Theta1(i,:),M11,M21);
        ddtheta21(i) = f2(Theta1(i,:),M11,M21);
        
        Theta2(i,:) = RK(Theta2(i-1,:),M12,M22,h);
        theta12(i) = Theta2(i,1);
        theta22(i) = Theta2(i,2);
        dtheta12(i) = Theta2(i,3);
        dtheta22(i) = Theta2(i,4);
        ddtheta12(i) = f1(Theta2(i,:),M12,M22);
        ddtheta22(i) = f2(Theta2(i,:),M12,M22);
    end
    
end

%time for problem 3
t3 = h*(0:1:m-1);

%solving problem 3
for i = 1:m
    if i == 1
        %initialize the matrix
        Theta3(1,:) = zeros(1,4);
    else
        %for each particular theta, 
        %carry out linear interplotion 
        M13 = interp1(theta1,M13d,Theta3(i-1,2),'linear');
        M23 = interp1(theta1,M23d,Theta3(i-1,2),'linear');
        
        Theta3(i,:) = RK(Theta3(i-1,:),M13,M23,h);
        theta13(i) = Theta3(i,1);
        theta23(i) = Theta3(i,2);
        dtheta13(i) = Theta3(i,3);
        dtheta23(i) = Theta3(i,4);
        ddtheta13(i) = f1(Theta3(i,:),M13,M23);
        ddtheta23(i) = f2(Theta3(i,:),M13,M23);
    end
    
end

%linear interplotation
theta = 0:pi/200:pi/2;
M13 = interp1(theta1,M13d,theta,'linear');
M23 = interp1(theta1,M23d,theta,'linear');

%draw figures
%problem 1
figure(1)
plot(t,theta11);
xlabel('t/s');
hold on
plot(t,theta21);
ylabel('\theta/rad');
title('\theta_1 & \theta_2, M_1 = 10, M_2 = 0');
legend('\theta_1','\theta_2');

figure(2)
plot(t,dtheta11);
xlabel('t/s');
hold on
plot(t,dtheta21);
ylabel('\omega/(rad/s)');
title('\omega_1 & \omega_2, M_1 = 10, M_2 = 0');
legend('\omega_1','\omega_2');

figure(3)
plot(t,ddtheta11);
xlabel('t/s');
hold on
plot(t,ddtheta21);
ylabel('\alpha/(rad/s^2)');
title('\alpha_1 & \alpha_2, M_1 = 10, M_2 = 0');
legend('\alpha_1','\alpha_2');

%problem 2
figure(4)
plot(t,theta12);
xlabel('t/s');
hold on
plot(t,theta22);
ylabel('\theta/rad');
title('\theta_1 & \theta_2, M_1 = 0, M_2 = 10');
legend('\theta_1','\theta_2');

figure(5)
plot(t,dtheta12);
xlabel('t/s');
hold on
plot(t,dtheta22);
ylabel('\omega/(rad/s)');
title('\omega_1 & \omega_2, M_1 = 0, M_2 = 10');
legend('\omega_1','\omega_2');

figure(6)
plot(t,ddtheta12);
xlabel('t');
hold on
plot(t,ddtheta22);
ylabel('\alpha/(rad/s^2)');
title('\alpha_1 & \alpha_2, M_1 = 0, M_2 = 10');
legend('\alpha_1','\alpha_2');

%problem 3
y1 = pi/2*ones(1,101);

x2 = 3*ones(1,101);
y2 = 0:pi/200:pi/2;

figure(7)
plot(t3,theta13,'g',t3,theta23,'b',t3,y1,'r-',x2,y2,'r-');
xlabel('t/s');
ylabel('\theta/rad');
title('\theta_1 & \theta_2, Synchronous Rotation');
legend('\theta_1','\theta_2');
axis([0 3.5 0 2]);

figure(8)
plot(t3,dtheta13);
xlabel('t/s');
ylabel('\omega/(rad/s)');
hold on
plot(t3,dtheta23);
title('\omega_1 & \omega_2, Synchronous Rotation');
legend('\omega_1','\omega_2');

figure(9)
plot(t3,ddtheta13);
xlabel('t/s');
ylabel('\alpha/(rad/s^2)');
hold on
plot(t3,ddtheta23);
title('\alpha_1 & \alpha_2, Synchronous Rotation');
legend('\alpha_1','\alpha_2');

%motor torque
figure(10)
y3 = zeros(1,m);
plot(theta,M13,theta,M23,theta,y3,'r');
xlabel('\theta/rad');
ylabel('Torque');
title('Motor Torque vs Rotational Angle');
legend('M_1','M_2');
axis([0 1.6 -30 35]);

%sub functions
%f1, angular acceleration for arm 1
function ddtheta1 = f1(Theta,M1,M2)

    D11 = 10 + 8*cos(Theta(2));
    D22 = 4;
    D12 = 4 + 4*cos(Theta(2));
    D21 = D12;
    D122 = -4*sin(Theta(2));
    D211 = -D122;
    D112 = -4*sin(Theta(2));
    
    C = D22*D11 - D12*D21;
    
    A0 = (D22*M1 - D12*M2)/C;
    A11 = D12*D211/C;
    A12 = -2*D22*D112/C;
    A22 = -D22*D122/C;
    
    ddtheta1 = A0 + A11*(Theta(3))^2 + ...
               A12*Theta(3)*Theta(4) + A22*(Theta(4))^2;

end

%f2, angular acceleration for arm 2
function ddtheta2 = f2(Theta,M1,M2)

    D11 = 10 + 8*cos(Theta(2));
    D22 = 4;
    D12 = 4 + 4*cos(Theta(2));
    D21 = D12;
    D122 = -4*sin(Theta(2));
    D211 = -D122;
    D112 = -4*sin(Theta(2));
    
    C = D22*D11 - D12*D21;
    
    B0 = (D11*M2 - D21*M1)/C;
    B11 = -D11*D211/C;
    B12 = 2*D21*D112/C;
    B22 = D21*D122/C;
    
    ddtheta2 = B0 + B11*(Theta(3))^2 + ...
               B12*Theta(3)*Theta(4) + B22*(Theta(4))^2;

end

%fstate, derivative of the 4 state variables
function dTheta = fstate(Theta,M1,M2)
    
    dTheta = zeros(1,4);
    dTheta(1) = Theta(3);
    dTheta(2) = Theta(4);
	dTheta(3) = f1(Theta,M1,M2);
	dTheta(4) = f2(Theta,M1,M2);
    
end

%RK, Runge-Kutta method
function Thetan = RK(Theta,M1,M2,h)
    
    k1 = fstate(Theta,M1,M2);
    k2 = fstate(Theta+h*k1/2,M1,M2);
    k3 = fstate(Theta+h*k2/2,M1,M2);
    k4 = fstate(Theta+h*k3,M1,M2);
    
    Thetan =  Theta + h*(k1+2*k2+2*k3+k4)/6;
    
end
