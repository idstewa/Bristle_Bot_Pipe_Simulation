%Vibration Robot Slip Right Equation Builder
%Ian Stewart
clear
clc
close all
%% Symbols
syms A a alpha_0 b c g I K k_cs L M mp mu_k R_r w varphi
syms theta(t) dtheta(t) ddtheta(t) x(t) dx(t) ddx(t) y(t) dy(t) ddy(t) dalpha(t) dbeta(t) t %#ok<NASGU> 
syms Theta DTheta DDTheta X DX DDX Y DY DDY DAlpha DDAlpha DBeta DDBeta dB ddB dd dD ddD Lambda1 Lambda2 Lambda3 Lambda4 slipRdir %#ok<NASGU> 
%% Slip Right Equations
%Constraints
T1 = sin(alpha_0)/cos(alpha_0);
T2 = 1/T1;
y = R_r*(T1*sin(theta) - cos(theta) + 1) + c*(sin(theta) + T1*cos(theta) - T1) + T1*x - T1*dB;
dy = subs(diff(y,t),[diff(theta,t) diff(x,t)],[dtheta dx]);
dalpha = (x + c*(cos(theta) - 1) + R_r*sin(theta) - dB)/(L*cos(alpha_0));
DDAlpha = subs(diff(dalpha,t),[theta diff(theta,t) x diff(x,t)],[Theta DTheta X DX]);
dbeta = (2*R_r*(cos(theta) - 1))/(L*sin(alpha_0)) - (R_r*sin(theta)+ c*(cos(theta) - 1) + x - dB)/(L*cos(alpha_0));
DDBeta = subs(diff(dbeta,t),[theta diff(theta,t) x diff(x,t)],[Theta DTheta X DX]);
%ddB is zero
%Lagrange
vp = a*w*sin(w*t);
dvp = diff(vp,t);
T = 0.5*M*(dx^2 + dy^2) + 0.5*mp*(b^2*dtheta^2 + vp^2 + 2*vp*(dx*sin(theta) + dy*cos(theta)) + 2*b*dtheta*(dx*sin(theta) + dy*cos(theta) + vp)) + 0.5*I*dtheta^2;
V = 0.5*K*(dalpha^2 + dbeta^2 + 2*theta^2 + 2*theta*(dalpha - dbeta) + 2*A*(dalpha + dbeta)) + M*g*(x*sin(varphi) + y*cos(varphi)) + 0.5*k_cs*theta^2;
La = T - V;
%Theta equation
dLadthetadot = diff(La,dtheta(t));
ddLadthetadotdt = diff(dLadthetadot,t);
dLadtheta = diff(La,theta(t));
La_theta = ddLadthetadotdt - dLadtheta == R_r*(-Lambda4)*cos(theta);
La_Theta = subs(La_theta,[theta diff(theta,t) dtheta diff(dtheta,t) x diff(x,t) dx diff(dx,t)],[Theta DTheta DTheta DDTheta X DX DX DDX]);
DDTheta = solve(La_Theta,DDTheta);
%X equation
dLadxdot = diff(La,dx(t));
ddLadxdotdt = diff(dLadxdot,t);
dLadx = diff(La,x(t));
La_x = ddLadxdotdt - dLadx == Lambda4;
La_X = subs(La_x,[theta diff(theta,t) dtheta diff(dtheta,t) x diff(x,t) dx diff(dx,t)],[Theta DTheta DTheta DDTheta X DX DX DDX]);
DDX = solve(La_X,DDX);
strDDX = replace(string(DDX),{'DTheta' 'Theta' 'DX' 'X'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)'}) %#ok<NOPTS> 
%Finish theta equation
DDTheta = simplify(subs(DDTheta));
strDDTheta = replace(string(DDTheta),{'DTheta' 'Theta' 'DX' 'X'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)'}) %#ok<NOPTS> 
%Lambdas
DDY = subs(diff(dy,t),[theta diff(theta,t) dtheta diff(dtheta,t) diff(dx,t)],[Theta DTheta DTheta DDTheta DDX]);
DAlpha = subs(dalpha,[theta x],[Theta X]);
DBeta = subs(dbeta,[theta x],[Theta X]);
c1 = L*sin(alpha_0);
c2 = L*cos(alpha_0);
f2 = M*DDY + mp*b*DDTheta*cos(Theta) - mp*b*DTheta^2*sin(Theta) - mp*vp*DTheta*sin(Theta) + mp*dvp*cos(Theta) + M*g*cos(varphi);
f3 = K*(DAlpha + Theta + A);
f4 = K*(DBeta - Theta + A);
Lambda3 = simplify(f4/(-c2*mu_k*slipRdir - c1));
Lambda4 = -abs(mu_k*Lambda3)*slipRdir;
f2 = subs(f2);
Lambda1 = (f2-Lambda3);
Lambda2 = (f3 - c1*Lambda1)/c2;
% strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta' 'DX' 'X' 't'},{'output(i,2)' 'output(i,1)' 'output(i,4)' 'output(i,3)' 'time(i,1)'}) %#ok<NASGU,NOPTS>
% strLambda2 = replace(string(Lambda2),{'DTheta' 'Theta' 'DX' 'X' 't'},{'output(i,2)' 'output(i,1)' 'output(i,4)' 'output(i,3)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
% strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta' 'DX' 'X' 't'},{'output(i,2)' 'output(i,1)' 'output(i,4)' 'output(i,3)' 'time(i,1)'}) %#ok<NASGU,NOPTS> 
% strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta' 'DX' 'X'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)'}) %#ok<NOPTS> 
% strLambda2 = replace(string(Lambda2),{'DTheta' 'Theta' 'DX' 'X'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)'}) %#ok<NOPTS> 
% strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta' 'DX' 'X'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)'}) %#ok<NOPTS> 
dD = x + c*(cos(theta) - 1) - R_r*sin(theta) - L*cos(alpha_0)*dbeta;
DD = subs(dD,[theta x],[Theta X]);
DDD = subs(diff(dD,t),[theta diff(theta,t) x diff(x,t)],[Theta DTheta X DX]);
strdD = replace(string(DD),{'Theta' 'X'},{'output(i,1)' 'output(i,3)'}) %#ok<NASGU,NOPTS> 
strddD = replace(string(DDD),{'DTheta' 'Theta' 'DX'},{'output(i,2)' 'output(i,1)' 'output(i,4)'}) %#ok<NOPTS,NASGU> 
strdD = replace(string(DD),{'Theta' 'X'},{'Z(1)' 'Z(3)'}) %#ok<NOPTS> 
strddD = replace(string(DDD),{'DTheta' 'Theta' 'DX'},{'Z(2)' 'Z(1)' 'Z(4)'}) %#ok<NOPTS> 

X = solve(DD==dd,X)
DX = solve(DDD==0,DX)
