%Author: JingXuan Yang
%Date: 2019.03.18
%Dynamics of Mechanical System HW1, part 1

%initial data 
phi1_d = [270 300 330 360 375 390 420 450 480 510 540 570 600];
Fp_d = [0 0.54 1.81 3.34 5.80 2.40 0.835 0.194 0.194 0.145 0.096 0.048 0];

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

%linear intepolation
phi1 = 1:1:720;
Fp = zeros(1,720);
Fp(270:600) = interp1(phi1_d,Fp_d,phi1(270:600),'linear');

%
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

%equivalent moment of inertia, Je
for i = 1:720
    
    %deal with abnormal conditions
    if phi1(i) == 90
        phi1(i) = 89.999;      
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
    
    %omega and angular acceleration
    dphi1(i) = 1;
    
    %phi2
    phi2(i) = asind(l1*sind(phi1(i))/l2);
    dphi2(i) = dphi1(i)*l1*cosd(phi1(i))/(l2*cosd(phi2(i)));
    ddphi2(i) = (ddphi1(i)*l1*cosd(phi1(i))-(dphi1(i))^2*l1*sind(phi1(i)) ...
                 +(dphi2(i))^2*l2*sind(phi2(i)))/(l2*cosd(phi2(i)));
    
    %omega2, alpha2         
    omega2(i) = -dphi2(i);
    alpha2(i) = ddphi2(i);
    
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
    
    %dJe/dphi1
    dJe(i) = 2*(m2*(S2dx(i)*S2ddx(i)+S2dy(i)*S2ddy(i))+m3*Cdy(i)*Cddy(i));
        
end %w.r.t for

%draw figures
figure(1)
plot(phi1,Fp);
xlabel('\phi_1/\circ');
ylabel('F_p/MPa');
title('Work Pressure vs Rotational Angle (linear intepolation)');

figure(2)
plot(phi1,Je);
yyaxis left
xlabel('\phi_1/\circ');
ylabel('J_e/kg\cdotm^2');
title('Equivalent Moment of Inertia and its Derivative vs Rotational Angle');
yyaxis right
plot(phi1,dJe);
ylabel('dJ_e/d\phi/(kg\cdotm^2/\circ)');
legend('J_e','dJ_e/d\phi');

figure(3)
plot(Bx,By);
hold on
plot(S2x,S2y,'r');
hold on
plot(Cx,Cy,'c');
axis equal

figure(4)
Bv = sqrt(Bdx.^2+Bdy.^2);
S2v = sqrt(S2dx.^2+S2dy.^2);
plot(phi1,Bv);
hold on
plot(phi1,S2v);

