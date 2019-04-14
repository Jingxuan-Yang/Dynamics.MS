%Author:  JingXuan Yang
%E-mail:  yangjingxuan@stu.hit.edu.cn
%Date:    2019.03.25
%Project: Dynamics of Mechanical System Homework 2
%Purpose: Dynamics of 2 DOF manipulater
%Note   : all angles in this script are in radius

% clear;  
% clc;

%initial data
%problem 1
M11 = 10;
M21 = 0;

%problem 2
M12 = 0;
M22 = 10;

% %problem 3
% M13d = [5, 30, 30, -4, -6, -25.05, -17.64,  0];
% M23d = [2.29, 13.88, 14.38, 0.4926, 2.628,  -10.28, -9.423, 0];
% rd = 180/pi;
% theta1 = [0, pi/80, 4*pi/80, pi/8,3*pi/8, 36*pi/80, 39*pi/80, pi/2];

% example in textbook
% M13d = [27.22,23.89,18.165,-6.58,-9.30,-11.39,-32.33,-23.63,-14.66];
% M23d = [12.58,12.99,14.195,2.094,3.10,3.797,-6.672,-7.34,-8.378];
% rd = 180/pi;
% dr = pi/180;
% theta1 = dr*[0, 15, 29.999, 30.001,45, 59.999,60.001,75,90];

%construct for debug
M13d = [5, 30, 30, -4, -6, -25.05, -17.64,  0];
M23d = [2.29, 13.88, 14.38, 0.4926, 2.628,  -10.28, -9.423, 0];
rd = 180/pi;
theta1 = [0, pi/80, 4*pi/80, pi/8,3*pi/8, 36*pi/80, 39*pi/80, pi/2];

%number of loop time
n = 1000;
m = 300;
%step
h = 0.01;

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

%calculating
%interplotation
theta = pi/180:pi/180:pi/2;
% M13 = interp1(theta1,M13d,theta,'linear');
% M23 = interp1(theta1,M23d,theta,'linear');
M13 = -12*sin(theta);
M23 = 4*sin(theta);

%time
t = 1:1:n;

for i = 1:n
    
    if i == 1
        theta11(i) = 0;
        theta21(i) = 0;
        dtheta11(i) = 0;
        dtheta21(i) = 0;
        
        theta12(i) = 0;
        theta22(i) = 0;
        dtheta12(i) = 0;
        dtheta22(i) = 0;        
    else
        [theta11(i),theta21(i),dtheta11(i),dtheta21(i)] = ...
            RK(theta11(i-1),theta21(i-1),dtheta11(i-1),dtheta21(i-1),M11,M21,h);
        ddtheta11(i) = f1(theta11(i),theta21(i),dtheta11(i),dtheta21(i),M11,M21);
        ddtheta21(i) = f2(theta11(i),theta21(i),dtheta11(i),dtheta21(i),M11,M21);
        
        [theta12(i),theta22(i),dtheta12(i),dtheta22(i)] = ...
            RK(theta12(i-1),theta22(i-1),dtheta12(i-1),dtheta22(i-1),M12,M22,h);
        ddtheta12(i) = f1(theta12(i),theta22(i),dtheta12(i),dtheta22(i),M12,M22);
        ddtheta22(i) = f2(theta12(i),theta22(i),dtheta12(i),dtheta22(i),M12,M22);        
    end
    
end

%time for problem 3
t3 = h*(1:1:m);
for i = 1:m
    if i == 1
        theta13(i) = 0;
        theta23(i) = 0;
        dtheta13(i) = 0;
        dtheta23(i) = 0;
    else
        [theta13(i),theta23(i),dtheta13(i),dtheta23(i)] = ...
            RK(theta13(i-1),theta23(i-1),dtheta13(i-1),dtheta23(i-1),M1(i),M2(i),h);
        ddtheta13(i) = f1(theta13(i),theta23(i),dtheta13(i),dtheta23(i),M1(i),M2(i));
        ddtheta23(i) = f2(theta13(i),theta23(i),dtheta13(i),dtheta23(i),M1(i),M2(i));
    end
    
end

%draw figure
%problem 1
figure(1)
plot(t,theta11,t,theta21);
xlabel('time');
ylabel('\theta_1');
title('\theta_1');

figure(2)
plot(t,theta21);
xlabel('time');
ylabel('\theta_2');
title('\theta_2');

figure(3)
plot(t,dtheta11);
xlabel('time');
ylabel('\omega_1');
title('\omega_1');

figure(4)
plot(t,dtheta21);
xlabel('time');
ylabel('\omega_2');
title('\omega_2');

figure(5)
plot(t,ddtheta11);
xlabel('time');
ylabel('\alpha_1');
title('\alpha_1');

figure(6)
plot(t,ddtheta21);
xlabel('time');
ylabel('\alpha_2');
title('\alpha_2');

%problem 2
figure(7)
plot(t,theta12);
xlabel('time');
ylabel('\theta_1');
title('\theta_1');

figure(8)
plot(t,theta22);
xlabel('time');
ylabel('\theta_2');
title('\theta_2');

figure(9)
plot(t,dtheta12);
xlabel('time');
ylabel('\omega_1');
title('\omega_1');

figure(10)
plot(t,dtheta22);
xlabel('time');
ylabel('\omega_2');
title('\omega_2');

figure(11)
plot(t,ddtheta12);
xlabel('time');
ylabel('\alpha_1');
title('\alpha_1');

figure(12)
plot(t,ddtheta22);
xlabel('time');
ylabel('\alpha_2');
title('\alpha_2');

%problem 3
figure(13)
plot(t3,theta13);
xlabel('time');
ylabel('\theta_1');
title('\theta_1');

figure(14)
plot(t3,theta23);
xlabel('time');
ylabel('\theta_2');
title('\theta_2');

figure(15)
plot(t3,dtheta13);
xlabel('time');
ylabel('\omega_1');
title('\omega_1');

figure(16)
plot(t3,dtheta23);
xlabel('time');
ylabel('\omega_2');
title('\omega_2');

figure(17)
plot(t3,ddtheta13);
xlabel('time');
ylabel('\alpha_1');
title('\alpha_1');

figure(18)
plot(t3,ddtheta23);
xlabel('time');
ylabel('\alpha_2');
title('\alpha_2');

%functions
function ddtheta1 = f1(~,theta2,dtheta1,dtheta2,M1,M2)

    D11 = 10 + 8*cos(theta2);
    D22 = 4;
    D12 = 4 + 4*cos(theta2);
    D21 = D12;
    D122 = -4*sin(theta2);
    D211 = -D122;
    D112 = -4*sin(theta2);
    
    C = D22*D11 - D12*D21;
    
    A0 = (D22*M1 - D12*M2)/C;
    A11 = D12*D211/C;
    A12 = -2*D22*D112/C;
    A22 = -D22*D122/C;
    
    ddtheta1 = A0 + A11*(dtheta1)^2 + A12*dtheta1*dtheta2 + A22*(dtheta2)^2;

end

function ddtheta2 = f2(~,theta2,dtheta1,dtheta2,M1,M2)

    D11 = 10 + 8*cos(theta2);
    D22 = 4;
    D12 = 4 + 4*cos(theta2);
    D21 = D12;
    D122 = -4*sin(theta2);
    D211 = -D122;
    D112 = -4*sin(theta2);
    
    C = D22*D11 - D12*D21;
    
    B0 = (D11*M2 - D21*M1)/C;
    B11 = -D11*D211/C;
    B12 = 2*D21*D112/C;
    B22 = D21*D122/C;
    
    ddtheta2 = B0 + B11*(dtheta1)^2 + B12*dtheta1*dtheta2 + B22*(dtheta2)^2;

end

function [th1n,th2n,dth1n,dth2n] = RK(theta1,theta2,dtheta1,dtheta2,M1,M2,h)

    c11 = h*f1(theta1,theta2,dtheta1,dtheta2,M1,M2);
    c21 = h*f2(theta1,theta2,dtheta1,dtheta2,M1,M2);
    
    d11 = h*dtheta1;
    d21 = h*dtheta2;
    
    c12 = h*f1(theta1+d11/2,theta2+d21/2,dtheta1+c11/2,dtheta2+c21/2,M1,M2);
    c22 = h*f2(theta1+d11/2,theta2+d21/2,dtheta1+c11/2,dtheta2+c21/2,M1,M2);
    
    d12 = h*(dtheta1+c11/2);
    d13 = h*(dtheta1+c12/2);
    d22 = h*(dtheta2+c21/2);
    d23 = h*(dtheta2+c22/2);
    
    c13 = h*f1(theta1+d12/2,theta2+d22/2,dtheta1+c12/2,dtheta2+c22/2,M1,M2);
    c23 = h*f2(theta1+d12/2,theta2+d22/2,dtheta1+c12/2,dtheta2+c22/2,M1,M2);

    d14 = h*(dtheta1+c13);
    d24 = h*(dtheta1+c23);
    
    c14 = h*f1(theta1+d13,theta2+d23,dtheta1+c13,dtheta2+c23,M1,M2);
    c24 = h*f2(theta1+d13,theta2+d23,dtheta1+c13,dtheta2+c23,M1,M2);
    
    th1n = theta1 + (d11 + 2*d12 + 2*d13 + d14)/6;
    th2n = theta2 + (d21 + 2*d22 + 2*d23 + d24)/6;
    
    dth1n = dtheta1 + (c11 + 2*c12 + 2*c13 + c14)/6;
    dth2n = dtheta2 + (c21 + 2*c22 + 2*c23 + c24)/6;
    
end


