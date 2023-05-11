%Vibration Robot Slip Both Equation Builder
%Ian Stewart
clear
clc
close all
%% Symbols
syms A a alpha_0 b c g I K k_cs L M mp mu_k mu_k_R mu_k_L R_r w varphi
syms theta(t) dtheta(t) ddtheta(t) x(t) dx(t) ddx(t) y(t) dy(t) ddy(t) dalpha(t) dbeta(t) t %#ok<NASGU> 
syms Theta DTheta DDTheta X DX DDX Y DY DDY DAlpha DBeta db dB ddb dd dD ddd Lambda1 Lambda2 Lambda3 Lambda4 slipRdir slipLdir %#ok<NASGU> 
%% Slip Both Equations
%Constraints
T1 = sin(alpha_0)/cos(alpha_0);
T2 = 1/T1;
dalpha = (y - c*sin(theta) + R_r*(cos(theta) - 1))/(L*sin(alpha_0));
DDAlpha = subs(diff(dalpha,t),[theta diff(theta,t) x diff(x,t)],[Theta DTheta X DX]);
dbeta = (c*sin(theta) + R_r*(cos(theta) - 1) - y)/(L*sin(alpha_0));
DDBeta = subs(diff(dbeta,t),[theta diff(theta,t) x diff(x,t)],[Theta DTheta X DX]);
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
La_theta = ddLadthetadotdt - dLadtheta == R_r*(Lambda2-Lambda4)*cos(theta);
La_Theta = subs(La_theta,[theta diff(theta,t) dtheta diff(dtheta,t) x diff(x,t) dx diff(dx,t) y diff(y,t) dy diff(dy,t)],[Theta DTheta DTheta DDTheta X DX DX DDX Y DY DY DDY]);
DDTheta = solve(La_Theta,DDTheta);
%X equation
dLadxdot = diff(La,dx(t));
ddLadxdotdt = diff(dLadxdot,t);
dLadx = diff(La,x(t));
La_x = ddLadxdotdt - dLadx == Lambda2+Lambda4;
La_X = subs(La_x,[theta diff(theta,t) dtheta diff(dtheta,t) x diff(x,t) dx diff(dx,t)],[Theta DTheta DTheta DDTheta X DX DX DDX]);
DDX = solve(La_X,DDX);
DDTheta = subs(DDTheta);
%Y equation
dLadydot = diff(La,dy(t));
ddLadydotdt = diff(dLadydot,t);
dLady = diff(La,y(t));
La_y = ddLadydotdt - dLady == 0;
La_Y = subs(La_y,[theta diff(theta,t) dtheta diff(dtheta,t) y diff(y,t) dy diff(dy,t)],[Theta DTheta DTheta DDTheta Y DY DY DDY]);
DDY = solve(La_Y,DDY);
strDDY = replace(string(DDY),{'DTheta' 'Theta' 'DX' 'X' 'DY' 'Y'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)' 'Z(6)' 'Z(5)'}) %#ok<NOPTS> 
%Finish x equation
DDX = subs(DDX);
strDDX = replace(string(DDX),{'DTheta' 'Theta' 'DX' 'X' 'DY' 'Y'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)' 'Z(6)' 'Z(5)'}) %#ok<NOPTS> 
%Finish theta equation
DDTheta = simplify(subs(DDTheta));
strDDTheta = replace(string(DDTheta),{'DTheta' 'Theta' 'DX' 'X' 'DY' 'Y'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(3)' 'Z(6)' 'Z(5)'}) %#ok<NOPTS> 
%Lambdas
DAlpha = subs(dalpha,[theta y],[Theta Y]);
DBeta = subs(dbeta,[theta y],[Theta Y]);
c1 = L*sin(alpha_0);
c2 = L*cos(alpha_0);
f3 = K*(DAlpha + Theta + A);
f4 = K*(DBeta - Theta + A);
Lambda1 = f3/(c1 + c2*mu_k_L*slipLdir);
Lambda2 = -abs(mu_k_L*Lambda1)*slipLdir;
Lambda3 = f4/(-c2*mu_k_R*slipRdir - c1);
Lambda4 = -abs(mu_k_R*Lambda3)*slipRdir;
% strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta' 'DY' 'Y' 't'},{'output(i,2)' 'output(i,1)' 'output(i,6)' 'output(i,5)' 'time(i,1)'})
% strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta' 'DY' 'Y' 't'},{'output(i,2)' 'output(i,1)' 'output(i,6)' 'output(i,5)' 'time(i,1)'})
% strLambda1 = replace(string(Lambda1),{'DTheta' 'Theta' 'DY' 'Y'},{'Z(2)' 'Z(1)' 'Z(6)' 'Z(5)'})
% strLambda3 = replace(string(Lambda3),{'DTheta' 'Theta' 'DY' 'Y'},{'Z(2)' 'Z(1)' 'Z(6)' 'Z(5)'})
dB = x + c*(cos(theta) - 1) + R_r*sin(theta) - L*cos(alpha_0)*dalpha;
DB = subs(dB,[theta x y],[Theta X Y]);
DDB = subs(diff(dB,t),[theta diff(theta,t) diff(x,t) diff(y,t)],[Theta DTheta DX DY]);
dD = x + c*(cos(theta) - 1) - R_r*sin(theta) - L*cos(alpha_0)*dbeta;
DD = subs(dD,[theta x y],[Theta X Y]);
DDD = subs(diff(dD,t),[theta diff(theta,t) diff(x,t) diff(y,t)],[Theta DTheta DX DY]);
strdB = replace(string(DB),{'Theta' 'X' 'Y'},{'output(i,1)' 'output(i,3)' 'output(i,5)'}) %#ok<NOPTS,NASGU> 
strddB = replace(string(DDB),{'DTheta' 'Theta' 'DX' 'DY'},{'output(i,2)' 'output(i,1)' 'output(i,4)' 'output(i,6)'}) %#ok<NOPTS,NASGU> 
strdD = replace(string(DD),{'Theta' 'X' 'Y'},{'output(i,1)' 'output(i,3)' 'output(i,5)'}) %#ok<NOPTS,NASGU> 
strddD = replace(string(DDD),{'DTheta' 'Theta' 'DX' 'DY'},{'output(i,2)' 'output(i,1)' 'output(i,4)' 'output(i,6)'}) %#ok<NOPTS,NASGU> 
strdB = replace(string(DB),{'Theta' 'X' 'Y'},{'Z(1)' 'Z(3)' 'Z(5)'}) %#ok<NOPTS> 
strddB = replace(string(DDB),{'DTheta' 'Theta' 'DX' 'DY'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(6)'}) %#ok<NOPTS> 
strdD = replace(string(DD),{'Theta' 'X' 'Y'},{'Z(1)' 'Z(3)' 'Z(5)'}) %#ok<NOPTS> 
strddD = replace(string(DDD),{'DTheta' 'Theta' 'DX' 'DY'},{'Z(2)' 'Z(1)' 'Z(4)' 'Z(6)'}) %#ok<NOPTS> 

Y1 = solve(DB==db,Y);
Y2 = solve(DD==dd,Y);
X = solve(Y1==Y2,X)
Y = subs(Y1)

DY1 = solve(DDB==ddb,DY);
DY2 = solve(DDD==0,DY);
DX = solve(DY1==DY2,DX)
DY = subs(DY1)
