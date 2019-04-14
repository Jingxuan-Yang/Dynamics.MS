%Author:  JingXuan Yang
%E-mail:  yangjingxuan@stu.hit.edu.cn
%Date:    2019.03.23
%Project: Dynamics of Mechanical System HW1
%Purpose: Solving the Motion State of Crank Slider

%initial data, 'd' means discrete
%in degree
phi1_d = [270 300 330 360 375 390 420 450 480 510 540 570 600];
%in Pa
Fp_d = 10^6*[0 0.54 1.81 3.34 5.80 2.40 0.835 0.194 0.194 0.145 0.096 0.048 0];

%length, in meter
l1 = 60*0.001;
l2 = 240*0.001;
ls2 = 80*0.001;
d = 100*0.001;

%moment of inertia, in kilogram times meter square
J1 = 0.49;
J2 = 0.049;

%gravity, in Newton
G1 = 206;
G2 = 19.6;
G3 = 9.8;
g = 9.8;

%mass, in kilogram
m1 = G1/g;
m2 = G2/g;
m3 = G3/g;

%coefficient of Mr w.r.t omega square
k = 0.001927;

%transform between degree and radius
rd = 180/pi;
dr = pi/180;

%step for Runge-Kutta method
h = 1*dr;

%in degree
phi1 = 1:1:720;
%in radius
phi1r = phi1*dr;

%linear interpolation of Fp
Fp = zeros(1,720);
Fp(270:600) = interp1(phi1_d,Fp_d,phi1(270:600),'linear');

%initialize the matrices
phi2 = zeros(1,720);
lks2 = zeros(1,720);
Je = zeros(1,720);
dphi1 = zeros(1,720);
ddphi1 = zeros(1,720);
dphi2 = zeros(1,720);
ddphi2 = zeros(1,720);
Bx = zeros(1,720);
By = zeros(1,720);
S2x = zeros(1,720);
S2y = zeros(1,720);
Cx = zeros(1,720);
Cy = zeros(1,720);
Bdx = zeros(1,720);
Bdy = zeros(1,720);
S2dx = zeros(1,720);
S2dy = zeros(1,720);
Cdx = zeros(1,720);
Cdy = zeros(1,720);
Bddx = zeros(1,720);
Bddy = zeros(1,720);
S2ddx = zeros(1,720);
S2ddy = zeros(1,720);
Cddx = zeros(1,720);
Cddy = zeros(1,720);
omega2 = zeros(1,720);
alpha2 = zeros(1,720);
dJe = zeros(1,720);
dJer = zeros(1,720);
Me = zeros(1,720);
Med = zeros(1,720);
Mer = zeros(1,720);

%calculating the Je and dJe/dphi, 
%all the angles in this 'for' loop is in degree
for i = 1:720
    
    %deal with abnormal conditions
    if phi1(i) == 90
        phi1(i) = 90.001;      
    end
    
    if phi1(i) == 270
        phi1(i) = 270.001;
    end
    
    if phi1(i) == 450
        phi1(i) = 450.001;
    end
    
    if phi1(i) == 630
        phi1(i) = 630.001;
    end
    
    %let omega = 1, since Je and dJe/dphi are not relevant to omega
    dphi1(i) = 1;
   
    %phi2, derived from phi1
    phi2(i) = asind(l1*sind(phi1(i))/l2);
    dphi2(i) = dphi1(i)*l1*cosd(phi1(i))/(l2*cosd(phi2(i)));
    ddphi2(i) = (ddphi1(i)*l1*cosd(phi1(i))-(dphi1(i))^2*l1*sind(phi1(i)) ...
                 +(dphi2(i))^2*l2*sind(phi2(i)))/(l2*cosd(phi2(i)));
    
    %omega2, alpha2         
    omega2(i) = -dphi2(i);
    alpha2(i) = -ddphi2(i);
    
    %length of KS2
    lks2(i) = sqrt(ls2^2+(l2*cosd(phi2(i))/cosd(phi1(i)))^2 ...
           - 2*ls2*l2*cosd(phi2(i))*cosd(phi2(i)+phi1(i))/cosd(phi1(i)));
    
    %equivalent moment of inertia, Je
    Je(i) = J1 + J2*(l1*cosd(phi1(i))/(l2*cosd(phi2(i))))^2 ... 
         + m2*(l1*cosd(phi1(i))*lks2(i)/(l2*cosd(phi2(i))))^2 ...
         + m3*(l1*sind(phi1(i))*(l1*cosd(phi1(i))+l2*cosd(phi2(i))) ...
         /(l2*cosd(phi2(i))))^2;
     
    %position
    Bx(i) = l1*sind(phi1(i));
    By(i) = l1*cosd(phi1(i));
    S2x(i) = l1*sind(phi1(i))-ls2*sind(phi2(i));
    S2y(i) = l1*cosd(phi1(i))+ls2*cosd(phi2(i));
    Cx(i) = 0;
    Cy(i) = l1*cosd(phi1(i))+l2*cosd(phi2(i));
    
    %velocity
    Bdx(i) = dphi1(i)*l1*cosd(phi1(i));
    Bdy(i) = -dphi1(i)*l1*sind(phi1(i));
    S2dx(i) = dphi1(i)*l1*cosd(phi1(i))-dphi2(i)*ls2*cosd(phi2(i));
    S2dy(i) = -dphi1(i)*l1*sind(phi1(i))-dphi2(i)*ls2*sind(phi2(i));
    Cdx(i) = 0;
    Cdy(i) = -dphi1(i)*l1*sind(phi1(i))-dphi2(i)*l2*sind(phi2(i));
    
    %acceleration
    Bddx(i) = ddphi1(i)*l1*cosd(phi1(i))-(dphi1(i))^2*l1*sind(phi1(i));
    Bddy(i) = -ddphi1(i)*l1*sind(phi1(i))-(dphi1(i))^2*l1*cosd(phi1(i));
    S2ddx(i) = ddphi1(i)*l1*cosd(phi1(i))-(dphi1(i))^2*l1*sind(phi1(i)) ...
            -ddphi2(i)*ls2*cosd(phi2(i))+(dphi2(i))^2*ls2*sind(phi2(i));
    S2ddy(i) = -ddphi1(i)*l1*sind(phi1(i))-(dphi1(i))^2*l1*cosd(phi1(i)) ...
            -ddphi2(i)*ls2*sind(phi2(i))-(dphi2(i))^2*ls2*cosd(phi2(i));
    Cddx(i) = 0;
    Cddy(i) = -ddphi1(i)*l1*sind(phi1(i))-(dphi1(i))^2*l1*cosd(phi1(i)) ...
            -ddphi2(i)*l2*sind(phi2(i))-(dphi2(i))^2*l2*cosd(phi2(i));
    
    %dJe/dphi1, in degree
    dJe(i) = 2*(J2*omega2(i)*alpha2(i)+ ...
             m2*(S2dx(i)*S2ddx(i)+S2dy(i)*S2ddy(i))+m3*Cdy(i)*Cddy(i));
    
    
end %w.r.t for

dphi1 = zeros(1,720);
ddphi1 = zeros(1,720);
dphi2 = zeros(1,720);
ddphi2 = zeros(1,720);

%calculating the real dynamics of this system
for i = 1:720
    
    if i == 1        
        %initial values
        dphi1(i) = 120; %rad/s
        ddphi1(i) = -k*(dphi1(i))^2/Je(1);
    else
        dphi1(i) = RungeKutta(phi1r(i-1),dphi1(i-1),Fp(i),Je(i),dJe(i),l1,l2,d,k,h);
        ddphi1(i) = dphi1(i)*frk(phi1r(i),dphi1(i),Fp(i),Je(i),dJe(i),l1,l2,d,k);
    end
    
    %deal with abnormal conditions
    if phi1(i) == 90
        phi1(i) = 90.001;
        phi1r(i) = phi1(i)*dr;
    end
    
    if phi1(i) == 270
        phi1(i) = 270.001;
        phi1r(i) = phi1(i)*dr;
    end
    
    if phi1(i) == 450
        phi1(i) = 450.001;
        phi1r(i) = phi1(i)*dr;
    end
    
    if phi1(i) == 630
        phi1(i) = 630.001;
        phi1r(i) = phi1(i)*dr;
    end
   
    %phi2, derived from phi1
    phi2(i) = asin(l1*sin(phi1r(i))/l2);
    dphi2(i) = dphi1(i)*l1*cos(phi1r(i))/(l2*cos(phi2(i)));
    ddphi2(i) = (ddphi1(i)*l1*cos(phi1r(i))-(dphi1(i))^2*l1*sin(phi1r(i)) ...
                 +(dphi2(i))^2*l2*sin(phi2(i)))/(l2*cos(phi2(i)));
    
    %omega2, alpha2         
    omega2(i) = -dphi2(i);
    alpha2(i) = -ddphi2(i);
         
    %position
    Bx(i) = l1*sin(phi1r(i));
    By(i) = l1*cos(phi1r(i));
    S2x(i) = l1*sin(phi1r(i))-ls2*sin(phi2(i));
    S2y(i) = l1*cos(phi1r(i))+ls2*cos(phi2(i));
    Cx(i) = 0;
    Cy(i) = l1*cos(phi1r(i))+l2*cos(phi2(i));
    
    %velocity
    Bdx(i) = dphi1(i)*l1*cos(phi1r(i));
    Bdy(i) = -dphi1(i)*l1*sin(phi1r(i));
    S2dx(i) = dphi1(i)*l1*cos(phi1r(i))-dphi2(i)*ls2*cos(phi2(i));
    S2dy(i) = -dphi1(i)*l1*sin(phi1r(i))-dphi2(i)*ls2*sin(phi2(i));
    Cdx(i) = 0;
    Cdy(i) = -dphi1(i)*l1*sin(phi1r(i))-dphi2(i)*l2*sin(phi2(i));
    
    %acceleration
    Bddx(i) = ddphi1(i)*l1*cos(phi1r(i))-(dphi1(i))^2*l1*sin(phi1r(i));
    Bddy(i) = -ddphi1(i)*l1*sin(phi1r(i))-(dphi1(i))^2*l1*cos(phi1r(i));
    S2ddx(i) = ddphi1(i)*l1*cos(phi1r(i))-(dphi1(i))^2*l1*sin(phi1r(i)) ...
            -ddphi2(i)*ls2*cos(phi2(i))+(dphi2(i))^2*ls2*sin(phi2(i));
    S2ddy(i) = -ddphi1(i)*l1*sin(phi1r(i))-(dphi1(i))^2*l1*cos(phi1r(i)) ...
            -ddphi2(i)*ls2*sin(phi2(i))-(dphi2(i))^2*ls2*cos(phi2(i));
    Cddx(i) = 0;
    Cddy(i) = -ddphi1(i)*l1*sin(phi1r(i))-(dphi1(i))^2*l1*cos(phi1r(i)) ...
            -ddphi2(i)*l2*sin(phi2(i))-(dphi2(i))^2*l2*cos(phi2(i));
    
    %driving torque    
    Med(i) =  (1/4)*pi*d^2*Fp(i)*l1*sin(phi1r(i)) ...
            *(l1*cos(phi1r(i))+l2*cos(phi2(i)))/(l2*cos(phi2(i)));
    
    %resistance torque    
    Mer(i) = k*(dphi1(i))^2;
    
    %equivalent torque, Me
    Me(i) = Med(i) - Mer(i);
        
end %w.r.t for

%draw figures
figure(1)
plot(phi1,Fp);
xlabel('\phi_1/\circ');
ylabel('F_p/Pa');
title('Work Pressure vs Rotational Angle (linear interpolation)');

figure(2)
plot(phi1,Je);
yyaxis left
xlabel('\phi_1/\circ');
ylabel('J_e/kg\cdotm^2');
title('Equivalent Moment of Inertia and its Derivative vs Rotational Angle');
yyaxis right
plot(phi1,dJe);
ylabel('dJ_e/d\phi/(kg\cdotm^2/rad)');
legend('J_e','dJ_e/d\phi');

figure(3)
plot(Bx,By);
hold on
plot(S2x,S2y,'r');
hold on
plot(Cx,Cy,'c');
axis equal
xlabel('x/m');
ylabel('y/m');
title('Position of Point B, S_2 and C');
legend('p_B','p_{S2}','p_C');

figure(4)
Bv = sqrt(Bdx.^2+Bdy.^2);
S2v = sqrt(S2dx.^2+S2dy.^2);
plot(phi1,Bv,'r');
yyaxis left
ylabel('v_B/(m/s)');
yyaxis right
plot(phi1,S2v,'b');
xlabel('\phi_1/\circ');
ylabel('v_{S2}/(m/s)');
title('Velocity of B & S_2 vs \phi_1');
legend('v_B','v_{S2}');

figure(5)
plot(phi1,phi2*rd);
xlabel('\phi_1/\circ');
ylabel('\phi_2/\circ');
title('Rotational Angle \phi_2 vs \phi_1');

figure(6)
plot(phi1,dphi1,'r');
yyaxis left
xlabel('\phi_1/\circ');
ylabel('\omega_1/(rad/s)');
yyaxis right
plot(phi1,-dphi2,'b');
ylabel('\omega_2/(rad/s)');
title('Rotational Angular velocity \omega_1 & \omega_2 vs \phi_1');
legend('\omega_1','\omega_2');

figure(7)
plot(phi1,ddphi1,'r');
hold on
plot(phi1,-ddphi2,'c');
xlabel('\phi_1/\circ');
ylabel('\alpha/(rad/s^2)');
title('Rotational Angular Acceleration \alpha_1 & \alpha_2 vs \phi_1');
legend('\alpha_1','\alpha_2');

figure(8)
plot(phi1,Cdy);
xlabel('\phi_1/\circ');
ylabel('v_C/(m/s)');
title('Velocity of the Slider');

figure(9)
plot(phi1,Cddy);
xlabel('\phi_1/\circ');
ylabel('a_c/(m/s^2)');
title('Acceleration of the Slider');

figure(10)
plot(phi1,Med);
hold on
plot(phi1,Mer,'r');
xlabel('\phi_1/\circ');
ylabel('Torque/(N\cdotm)');
title('Driving and Resistance Torque vs Rotational Angle');
legend('M_d','M_r');

figure(11)
plot(phi1,Cy);
xlabel('\phi_1/\circ');
ylabel('Cy/m');
title('Position of Point C (Slider)');

%frk, function for runge kutta calculation
%dphin means dphi next
function dphi1n = frk(phi1r,dphi1,Fp,Je,dJe,l1,l2,d,k)
    
    phi2 = asin(l1*sin(phi1r)/l2);
    Me = (1/4)*pi*d^2*Fp*l1*sin(phi1r) ...
            *(l1*cos(phi1r)+l2*cos(phi2))/(l2*cos(phi2)) ...
            - k*(dphi1)^2;
    
    dphi1n = (Me-0.5*dJe*dphi1^2)/(Je*dphi1);

end

%function for runge-kutta method
%omegan means omega next
function omegan = RungeKutta(phi1r,dphi1,Fp,Je,dJe,l1,l2,d,k,h)

    k1 = frk(phi1r,dphi1,Fp,Je,dJe,l1,l2,d,k);
    k2 = frk(phi1r+h/2,dphi1+k1*h/2,Fp,Je,dJe,l1,l2,d,k);
    k3 = frk(phi1r+h/2,dphi1+k2*h/2,Fp,Je,dJe,l1,l2,d,k);
    k4 = frk(phi1r+h,dphi1+k3*h,Fp,Je,dJe,l1,l2,d,k);
    
    omegan = dphi1 + h*(k1+2*k2+2*k3+k4)/6;

end

