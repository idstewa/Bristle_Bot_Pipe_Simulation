%New Sim
%Ian Stewart
clear
clc
close all
%% Main
SLIP = zeros(1,10);
TIME = zeros(1,10);
VELOCITY = zeros(1,10);
pos = 0;
rpm = [4215.95 4908.81 5255.24 5601.67 6294.53 6987.39 7680.25 8373.11 9065.97 9758.83];
for j = 1:1%0
    RPM = 6000%rpm(j)
    FullSim = Slip(RPM);
    pos = pos + 1;
    SLIP(1,pos) = FullSim(2,end)-FullSim(2,1);
    TIME(1,pos) = FullSim(1,end)-FullSim(1,1);
    VELOCITY(1,pos) = -SLIP(1,pos)/TIME(1,pos);
end
velocity = [-0.000783333 0 -0.000675 0.000191667 0.075 0.101 0.142 0.143 0.158 0.148];
plot(rpm,velocity,'-ob',rpm,VELOCITY,'--^r','LineWidth',2,'MarkerSize',10)
xlabel('Motor Frequency [RPM]','FontSize',12,'interpreter','latex')
ylabel('Avg Velocity ($\dot{x}$) [m/sec]','FontSize',12,'interpreter','latex')
legend({'Experimental Responce' 'Simulated Responce'},'FontSize',12,'interpreter','latex')
grid on

warning('off','all')
figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(2,:))
title('Axial Travel (x)')
ylabel('Linear Displacement [m]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(3,:));
title('Axial Velocity (dx)')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(4,:))
title('Radial Travel (y)')
ylabel('Linear Displacement [m]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(5,:));
title('Radial Velocity (dx)')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(6,:))
title('Theta1')
ylabel('Anglular Displacement [rad]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(7,:));
title('DTheta1')
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/s]')

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(8,:))
title('Theta2')
ylabel('Anglular Displacement [rad]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(9,:));
title('DTheta2')
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/s]')

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(18,:))
title('Travel of Theta2 Bristle Tip (B)')
ylabel('Linear Displacement [m]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(19,:))
title('Velocity of Theta2 Bristle Tip (dB)')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(20,:))
title('Travel of Theta1 Bristle Tip (D)')
ylabel('Linear Displacement [m]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(21,:))
title('Velocity of Theta1 Bristle Tip (dD)')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

figure
plot(FullSim(1,:),FullSim(10,:),FullSim(1,:),FullSim(14,:))
title('Normal Forces')
xlabel('Time [sec]')
ylabel('Friction Force [N]')
legend({'Normal Force Theta1' 'Normal Force Theta2'})

[a,g,K,L,M,mp,mu,R,Rb,Rp,T,theta0,w] = params(RPM);

figure
subplot(2,1,1)
plot(FullSim(1,:),FullSim(11,:),FullSim(1,:),FullSim(12,:),FullSim(1,:),FullSim(13,:))
title('Theta1 Friction Force')
ylabel('Friction Force [N]')
subplot(2,1,2)
plot(FullSim(1,:),FullSim(15,:),FullSim(1,:),FullSim(16,:),FullSim(1,:),FullSim(17,:))
title('Theta2 Friction Force')
xlabel('Time [sec]')
ylabel('Friction Force [N]')
%% Slip
function FullSim = Slip(RPM)
    thetar = acos(((28.3/2)-(20.43/2))/7.06);
    IC = [thetar 0 0 0];

    [a,g,K,L,M,mp,mu_big,mu_small,R,Rb,Rp,T,theta0,w] = params(RPM); %#ok<ASGLU> 
    tspan = 0:T/1000:1;
    [time,output] = ode45(@(t,Z) SlipODE(t,Z,a,g,K,L,M,mp,mu_big,mu_small,R,theta0,w),tspan,IC);
    
    y = zeros(1,length(time));
    dy = zeros(1,length(time));
    theta2 = zeros(1,length(time));
    dtheta2 = zeros(1,length(time));
    N1 = zeros(1,length(time));
    N2 = zeros(1,length(time));
    f1 = zeros(1,length(time));
    lambda1 = zeros(1,length(time));
    lambda3 = zeros(1,length(time));
    f2 = zeros(1,length(time));
    lambda2 = zeros(1,length(time));
    lambda4 = zeros(1,length(time));
    db = zeros(1,length(time));
    ddb = zeros(1,length(time));
    dd = zeros(1,length(time));
    ddd = zeros(1,length(time));

    for i = 1:length(time)
        theta1 = output(i,1); dtheta1 = output(i,2); x = output(i,3); dx = output(i,4);
        y(1,i) = L*cos(theta1) + Rb;
        dy(1,i) = -L*dtheta1*sin(theta1);
        theta2(1,i) = acos(R-cos(theta1));
        dtheta2(1,i) = -dtheta1*sin(theta1)/sin(theta2(1,i));
        if sign(dx-L*dtheta1*cos(theta1))>=0
            mu1 = mu_small;
        else
            mu1 = mu_big;
        end
        if sign(dx-L*dtheta2(1,i)*cos(theta2(1,i)))>=0
            mu2 = mu_small;
        else
            mu2 = mu_big;
        end
        N1(1,i) = K*(theta1-theta0)/(L*sin(theta1) + L*cos(theta1)*mu1*sign(dx-L*dtheta1*cos(theta1)));
        N2(1,i) = K*(theta2(1,i)-theta0)/(L*cos(theta2(1,i))*mu2*sign(dx-L*dtheta2(1,i)*cos(theta2(1,i)))-L*sin(theta2(1,i)));
        f1(1,i) = -mu1*N1(1,i)*sign(dx-L*dtheta1*cos(theta1));
        lambda3(1,i) = abs((N1(1,i)*L*sin(theta1)-K*(theta1-theta0))/(L*cos(theta1)))*sign(f1(1,i));
        f2(1,i) = -mu2*N2(1,i)*sign(dx-L*dtheta2(1,i)*cos(theta2(1,i)));
        lambda4(1,i) = abs((N2(1,i)*L*sin(theta2) + K*(theta2 - theta0))/(L*cos(theta2)))*sign(f2(1,i));
        db(1,i) = x-L*sin(theta2(1,i));
        ddb(1,i) = dx-L*dtheta2(1,i)*cos(theta2(1,i));
        dd(1,i) = x-L*sin(theta1);
        ddd(1,i) = dx-L*dtheta1*cos(theta1);
    end

    FullSim = [time';output(:,3:4)';y;dy;output(:,1:2)';theta2;dtheta2;N1;f1;lambda1;lambda3;N2;f2;lambda2;lambda4;db;ddb;dd;ddd];
end
%% ODE Slip
function DZDT = SlipODE(t,Z,a,g,K,L,M,mp,mu_big,mu_small,R,theta0,w)
    theta1 = Z(1); dtheta1 = Z(2); dx = Z(4);
    theta2 = acos(R-cos(theta1));
    dtheta2 = -dtheta1*sin(theta1)/sin(theta2);
    if sign(dx-L*dtheta1*cos(theta1))>=0
        mu1 = mu_small;
    else
        mu1 = mu_big;
    end
    if sign(dx-L*dtheta2*cos(theta2))>=0
        mu2 = mu_small;
    else
        mu2 = mu_big;
    end
    N1 = K*(theta1-theta0)/(L*sin(theta1)+L*cos(theta1)*mu1*sign(dx-L*dtheta1*cos(theta1)));
    N2 = K*(theta2-theta0)/(L*cos(theta2)*mu2*sign(dx-L*dtheta2*cos(theta2))-L*sin(theta2));
    DZDT = zeros(4,1);
    %Theta1 dot
    DZDT(1) = Z(2);
    %Theta1 double dot
    DZDT(2) = -(N1+N2-M*g-mp*a*w^2*cos(w*t)-M*L*dtheta1^2*cos(theta1))/(M*L*sin(theta1));
    %X dot
    DZDT(3) = Z(4);
    %X double dot
    DZDT(4) = (-mu1*N1*sign(dx-L*dtheta1*cos(theta1))-mu2*N2*sign(dx-L*dtheta2*cos(theta2)))/M;
end
%% Parameters
function [a,g,K,L,M,mp,mu_big,mu_small,R,Rb,Rp,T,theta0,w] = params(RPM)
    a = 1.85/1000;
    g = 9.81;
    K = 2*(0.0246);
    L = 7.06/1000;
    M = 11.218/1000;
    mp = 0.942/1000;
    mu_big = 0.19;
    mu_small = 0.04;
    Rb = (20.43/2)/1000;
    Rp = (28.3/2)/1000;
    R = 2*(Rp-Rb)/L;
    theta0 = pi/4;
    Hz = RPM/60;
    T = 1/Hz;
    w = 2*pi*Hz;
end