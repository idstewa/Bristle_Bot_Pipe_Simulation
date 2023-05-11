%Vibration Robot Stick Mode Equation Builder
%Ian Stewart
clear
clc
close all
%% Symbols
syms A a a1 a2 alpha_0 b c CG g I K k_cs L M mm mp R_r w varphi
syms ddalpha(t) ddbeta(t) theta(t) dtheta(t) ddtheta(t) x(t) dx(t) ddx(t) y(t) dy(t) ddy(t) dalpha(t) dbeta(t) t %#ok<NASGU> 
syms Theta DTheta DDTheta X DX DDX Y DY DDY DAlpha DBeta dB dD Lambda1 Lambda2 Lambda3 Lambda4 %#ok<NASGU> 
%% Stick Equations
%Constraints
T1 = sin(alpha_0)/cos(alpha_0);
T2 = 1/T1;
x = (R_r*T2 - c)*(cos(theta) - 1) + 0.5*(dB + dD);
dx = subs(diff(x,t),diff(theta,t),dtheta);
y = (R_r*T1 + c)*sin(theta)+ 0.5*T1*(dD - dB);
dy = subs(diff(y,t),diff(theta,t),dtheta);
dalpha = (R_r/L)*(((cos(theta) - 1)/sin(alpha_0)) + (sin(theta)/cos(alpha_0))) + (dB - dD)/(2*L*cos(alpha_0));
ddalpha = subs(diff(dalpha,t),[theta diff(theta,t)],[Theta DTheta]);
dbeta = (R_r/L)*(((cos(theta) - 1)/sin(alpha_0)) - (sin(theta)/cos(alpha_0)))+ (dD - dB)/(2*L*cos(alpha_0));
ddbeta = subs(diff(dbeta,t),[theta diff(theta,t)],[Theta DTheta]);
%ddB and ddD are zero
%Lagrange
vp = a*w*sin(w*t);
dvp =diff(vp,t);
T = 0.5*M*(dx^2 + dy^2) + 0.5*mp*(b^2*dtheta^2 + vp^2 + 2*vp*(dx*sin(theta) + dy*cos(theta)) + 2*b*dtheta*(dx*sin(theta) + dy*cos(theta) + vp)) + 0.5*I*dtheta^2;
V = 0.5*K*(dalpha^2 + dbeta^2 + 2*theta^2 + 2*theta*(dalpha - dbeta) + 2*A*(dalpha + dbeta)) + M*g*(x*sin(varphi) + y*cos(varphi)) + 0.5*k_cs*theta^2;
La = T - V;
%Theta equation
dLadthetadot = diff(La,dtheta(t));
ddLadthetadotdt = diff(dLadthetadot,t);
dLadtheta = diff(La,theta(t));
La_theta = ddLadthetadotdt - dLadtheta == 0;
La_Theta = subs(La_theta,[theta diff(theta,t) dtheta diff(dtheta,t)],[Theta DTheta DTheta DDTheta]);
DDTheta = solve(La_Theta,DDTheta);
strDDTheta = replace(string(DDTheta),{'DTheta' 'Theta'},{'Z(2)' 'Z(1)'}) %#ok<NOPTS> 
%Lambdas\
DDX = subs(diff(dx,t),[theta diff(theta,t) dtheta diff(dtheta,t)],[Theta DTheta DTheta DDTheta]);
DDY = subs(diff(dy,t),[theta diff(theta,t) dtheta diff(dtheta,t)],[Theta DTheta DTheta DDTheta]);
DAlpha = subs(dalpha,theta,Theta);
DBeta = subs(dbeta,theta,Theta);
c1 = L*sin(alpha_0);
c2 = L*cos(alpha_0);
f1 = M*DDX + mp*b*DDTheta*sin(Theta) + mp*b*DTheta^2*cos(Theta) + mp*vp*DTheta*cos(Theta) + mp*dvp*sin(Theta) + M*g*sin(varphi);
f2 = M*DDY + mp*b*DDTheta*cos(Theta) - mp*b*DTheta^2*sin(Theta) - mp*vp*DTheta*sin(Theta) + mp*dvp*cos(Theta) + M*g*cos(varphi);
f3 = K*(DAlpha + Theta + A);
f4 = K*(DBeta - Theta + A);
Lambda4 = (c2*f1 + c1*f2 - f3 + f4) / (2*c2);
Lambda2 = f1 - Lambda4;
Lambda3 = (c2*Lambda4 - f4) / c1;
Lambda1 = f2 - Lambda3;
strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta' 't'},{'output(i,2)' 'output(i,1)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
strLambda2 = replace(string(Lambda2),{'DTheta' 'Theta' 't'},{'output(i,2)' 'output(i,1)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta' 't'},{'output(i,2)' 'output(i,1)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
strLambda4 = replace(string(Lambda4),{'DTheta' 'Theta' 't'},{'output(i,2)' 'output(i,1)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta'},{'Z(2)' 'Z(1)'}) %#ok<NOPTS> 
strLambda2 = replace(string(Lambda2),{'DTheta' 'Theta'},{'Z(2)' 'Z(1)'}) %#ok<NOPTS> 
strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta'},{'Z(2)' 'Z(1)'}) %#ok<NOPTS> 
strLambda4 = replace(string(Lambda4),{'DTheta' 'Theta'},{'Z(2)' 'Z(1)'}) %#ok<NOPTS> 






