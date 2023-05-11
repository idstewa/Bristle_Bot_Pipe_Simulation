%Lateral Slip Equation Builder
%Ian STewart
clear
clc
%% Symbols
syms a alpha0 AD BC CG g I k L M mp mu phi w
syms dx(t) y(t) dy(t) theta(t) dtheta(t) alpha(t) dalpha(t) beta(t) dbeta(t) t
syms DX DDX DY DDY Theta DTheta DDTheta Alpha DAlpha DDAlpha Beta DBeta DDBeta Lambda4
%% Stick Equations
%Lagrange
T = 0.5*M*(dx^2 + dy^2) - mp*a*w*(dx*cos(w*t + phi) - dy*sin(w*t + phi)) + 0.5*I*dtheta^2;
V = M*g*y + 0.5*k*(alpha^2 + beta^2 + 2*theta^2 + 2*alpha*(theta - alpha0) - 2*beta*(theta + alpha0));
La = T - V;
%X equation
dLadxdot = diff(La,dx(t));
ddLadxdotdt = diff(dLadxdot,t);
La_x = ddLadxdotdt;
f1 = subs(La_x,diff(dx,t),DDX);
x = 0.5*AD - (L*sin(alpha) + 0.5*BC*cos(theta) - CG*sin(theta));
dx = subs(diff(x,t),[diff(theta,t) diff(alpha,t)],[dtheta dalpha]);
DDX = subs(diff(dx,t),[theta diff(theta,t) dtheta diff(dtheta,t) alpha diff(alpha,t) dalpha diff(dalpha,t)],[Theta DTheta DTheta DDTheta Alpha DAlpha DAlpha DDAlpha]);
f1 = subs(f1);
%Y equation
dLadydot = diff(La,dy(t));
ddLadydotdt = diff(dLadydot,t);
dLady = diff(La,y(t));
La_y = ddLadydotdt - dLady;
f2 = subs(La_y,diff(dy,t),DDY);
y = L*cos(alpha) + 0.5*BC*sin(theta) + CG*cos(theta);
dy = subs(diff(y,t),[diff(theta,t) diff(alpha,t)],[dtheta dalpha]);
DDY = subs(diff(dy,t),[theta diff(theta,t) dtheta diff(dtheta,t) alpha diff(alpha,t) dalpha diff(dalpha,t)],[Theta DTheta DTheta DDTheta Alpha DAlpha DAlpha DDAlpha]);
f2 = subs(f2);
%Theta equation
dLadthetadot = diff(La,dtheta(t));
ddLadthetadotdt = diff(dLadthetadot,t);
dLadtheta = diff(La,theta(t));
La_theta = ddLadthetadotdt - dLadtheta;
f3 = subs(La_theta,[theta diff(dtheta,t) alpha beta],[Theta DDTheta Alpha Beta]);
dtheta = L*(cos(alpha)*dalpha + cos(beta)*dbeta)/(BC*sin(theta));
DDTheta = subs(diff(dtheta,t),[theta diff(theta,t) alpha diff(alpha,t) dalpha diff(dalpha,t) beta diff(beta,t) dbeta diff(dbeta,t)],[Theta DTheta Alpha DAlpha DAlpha DDAlpha Beta DBeta DBeta DDBeta]);
f1 = subs(f1);
f2 = subs(f2);
f3 = subs(f3);
%Beta equation
dLadbeta = diff(La,beta(t));
La_beta = -dLadbeta;
f5 = subs(La_beta,[theta beta],[Theta Beta]);
dbeta = dalpha*(sin(alpha)*sin(theta) - cos(alpha)*cos(theta))/(sin(beta)*sin(theta) + cos(beta)*cos(theta));
DDBeta = subs(diff(dbeta,t),[theta diff(theta,t) alpha diff(alpha,t) dalpha diff(dalpha,t) beta diff(beta,t)],[Theta DTheta Alpha DAlpha DAlpha DDAlpha Beta DBeta]);
f1 = subs(f1);
f2 = subs(f2);
f3 = subs(f3);
%Lambdas
c1 = -(0.5*BC*cos(Theta) - CG*sin(Theta));
c2 = -(0.5*BC*cos(Theta) + CG*sin(Theta));
c3 = 0.5*BC*cos(Theta) + CG*sin(Theta);
c4 = 0.5*BC*cos(Theta) - CG*sin(Theta);
c5 = L*sin(Beta);
c6 = -L*cos(Beta);
Lambda3 = (f5 - c6*Lambda4)/c5;
Lambda2 = f1 - Lambda4;
Lambda1 = f2 - Lambda3;
F3 = c1*Lambda1 + c2*Lambda2 + c3*Lambda3 + c4*Lambda4;
Lambda4 = solve(f3==F3,Lambda4);
Lambda3 = (f5 - c6*Lambda4)/c5;
Lambda2 = f1 - Lambda4;
Lambda1 = f2 - Lambda3;
f4 = k*(Alpha + Theta - alpha0) == L*sin(Alpha)*Lambda1 + L*cos(Alpha)*Lambda2;
DDAlpha = solve(f4,DDAlpha);
Lambda4 = subs(Lambda4);
Lambda3 = subs(Lambda3);
Lambda2 = subs(Lambda2);
Lambda1 = subs(Lambda1);
