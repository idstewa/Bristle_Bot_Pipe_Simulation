%Ian Stewart
%Vibration Robot Plots
clear
clc
close all
%% Voltage Vs RPM
V = [0.5 1.0 1.2 1.4 1.5 1.6 1.8 2.2 2.4 2.6 2.8 3.0];
RPM = [1495 3475 4245 4986 5357 5732 6383 7826 8491 9137 9730 10056];
x = 0.5:0.1:3;
y = 3464.3*x+58.279;
figure
plot(V,RPM,'ob',x,y,'-r','LineWidth',2,'MarkerSize',10)
grid on
caption = sprintf('RPM = 3464.3V+58.279');
text(0.6,9500,caption,'FontSize',16,'FontWeight','bold');
xlabel('Voltage (V) [V]')
ylabel('RPM')
%% Voltage Vs Speed
%Full Plot
V = [0.5 1 1.2 1.4 1.5 1.6 1.8 2 2.2 2.4 2.6 2.8 3];
avgxdot1 = [0 -0.000158333 -0.000783333 0 -0.00068 0.00019 0.075 0.101 0.142 0.143 0.158 0.148 0.16];
RPM = [1790.429 3522.6 4215.4 4908.3 5254.7 5601.2 6294.0 6986.9 7679.7 8372.6 9065.5 9758.3 10451.179];
avgxdot2 = [0 -0.0002 -0.0006 0.008 0.006 0.023 0.027 0.042 0.062 0.050 0.050 0.12 0.11];
figure
t = tiledlayout(1,1);
ax2 = axes(t);
% plot(ax2,RPM,1,'w')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.YColor = 'none';
axis(ax2,[1790.429 10451.179 -0.02 0.18])
ax1 = axes(t);
plot(ax1,V,avgxdot1,'-ob','LineWidth',2,'MarkerSize',10)
ax1.Box = 'off';
ax2.Box = 'off';
hold on
plot(ax1,V,avgxdot2,'--^r','LineWidth',2,'MarkerSize',10)
grid on
xlabel(ax1,'Voltage (V) [V]','FontSize',12,'interpreter','latex')
xlabel(ax2,'Motor Frequency ($\omega$) [RPM]','FontSize',12,'interpreter','latex')
ylabel(ax1,'Avg Velocity ($\dot{x}$) [m/sec]','FontSize',12,'interpreter','latex')
legend({'Larger Diameter Pipe' 'Smaller Diameter Pipe'},'FontSize',12,'interpreter','latex')
%Low Voltage Colse Up
V = [0.5 1 1.2 1.4 1.5 1.6];
avgxdot1 = [0 -0.000158333 -0.000783333 0 -0.00068 0.00019];
figure
t = tiledlayout(1,1);
ax4 = axes(t);
% plot(ax2,RPM,-0.01,'w')
ax4.XAxisLocation = 'top';
ax4.YAxisLocation = 'right';
ax4.YColor = 'none';
axis(ax4,[1443.999 5601.2 -0.02 0])
ax3 = axes(t);
plot(ax3,V,avgxdot1,'-ob','LineWidth',2,'MarkerSize',10)
grid on
xlabel(ax3,'Voltage (V) [V]','FontSize',12,'interpreter','latex')
xlabel(ax4,'Motor Frequency ($\omega$) [RPM]','FontSize',12,'interpreter','latex')
ylabel(ax3,'Avg Velocity ($\dot{x}$) [m/sec]','FontSize',12,'interpreter','latex')
ax3.YAxis.Exponent = 0;
%% Stiffness
angle = [0.8369 0.8659 0.8691 0.8921 0.8865 0.9017 0.9179 0.9439];
FL = [0.006819 0.007682 0.007558 0.008190 0.008089 0.008373 0.009241 0.009275];
err = [0.00015 0.00015 0.00015 0.00015 0.00015 0.00015 0.00015 0.00015];
x = 0.8369:0.0001:0.9439;
y = 0.0246*x-0.0137;
figure
errorbar(angle,FL,err,'o')
hold on
plot(x,y,'-r','LineWidth',2)
grid on
caption = sprintf('T=0.0246$\phi$-0.0137');
text(0.83,9.25*10^-3,caption,'FontSize',16,'FontWeight','bold','interpreter','latex');
xlabel('Cilia-Substrate Angle ($\phi$) [rad]','FontSize',12,'interpreter','latex')
ylabel('Torque (T) [Nm]','FontSize',12,'interpreter','latex')
% %% Friction
% XcellData = xlsread('FrictionData.xlsx'); %#ok<XLSRD> 
% figure
% plot(XcellData(:,1),XcellData(:,2),'LineWidth',2)
% grid on
% xlim([144 156])
% xlabel('Time (t) [sec]')
% ylabel('Force (F) [N]')
%% Friction New
XcellData = xlsread('Friction_Test_Data.xlsx','UpSideDown1'); %#ok<XLSRD>
figure
subplot(2,1,1)
plot(XcellData(2:1002,1),XcellData(2:1002,2),'-b','LineWidth',2)
grid on
xlim([0 10])

ylabel('Force (F) [N]','FontSize',12,'interpreter','latex')
title('Bristles Leading','FontSize',12,'interpreter','latex')
XcellData = xlsread('Friction_Test_Data.xlsx','RightSideUp1'); %#ok<XLSRD>
subplot(2,1,2)
plot(XcellData(2:1002,1),XcellData(2:1002,2),'-b','LineWidth',2)
grid on
xlim([0 10])
title('Bristles Trailing','FontSize',12,'interpreter','latex')
xlabel('Time (t) [sec]','FontSize',12,'interpreter','latex')
ylabel('Force (F) [N]','FontSize',12,'interpreter','latex')
%% Sim vs Exp
% XcellData = xlsread('1_Vibration Robot Data.xlsx','Sim2'); %#ok<XLSRD>
% figure
% plot(XcellData(1:13,1),XcellData(1:13,2),'-ob',XcellData(1:46,3),XcellData(1:46,12),'--^r','LineWidth',2,'MarkerSize',5)
% title('Experiment vs Simulation')
% xlabel('Motor Frequency [RPM]')
% ylabel('Avg Velocity ($\dot{x}$) [m/sec]')
% legend({'Experimental Responce' 'Simulated Responce'})
% grid on
%% Lateral Slip
volts1 = [12 6 2 1.7 1.5 1];
Velocity1 = [-3.575298217 -1.621048841 -0.317339436 0.104999182 0.056069161 0];
volts2 = [12 6 2.5 1.7 1.5 1];
Velocity2 = [3.636102608 1.118800803 0.277590462 -0.161901355 -0.153701698 0];
figure
plot(volts1,Velocity1,'-ob',volts2,Velocity2,'--^r','LineWidth',2,'MarkerSize',5)
xlabel('Voltage (V) [V]','FontSize',12,'interpreter','latex')
ylabel('Avg Velocity ($\dot{s}$) [m/sec]','FontSize',12,'interpreter','latex')
legend({'Clockwise' 'Counter Clockwise'},'FontSize',12,'interpreter','latex')
grid on
