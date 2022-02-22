syms L r rp h b R
%Optimization
d = 7850;
t = transpose([L r rp h b R]) %t=transpose([r lb rj rp])
fun = @(t)((pi*t(2)^2*t(4)*d)-(pi*t(3)^2*t(4)*d) + (6*t(1)*t(5)*t(4)*d))
x0=[0.020;0.015;0.002;0.010;0.010; 0.050]
A = [-1.000 0.000 0.000 0.000 0.000 0.000];
B = [-0.035] ; 
Aeq=[]; 
Beq=[];
 
% L = 0.02;
% r = 0.015;
% rp = 0.002;
% h = 0.010;  
% b = 0.010;

xLb=[0.01; 0.005; 0.001; 0.01; 0.01; 0.03];
xUb=[0.05; 0.025; 0.003; 0.03; 0.03; 0.05];

Nonlincon = @nonlincon_1;

x=fmincon(fun,x0,A,B,Aeq,Beq,xLb,xUb,Nonlincon)

disp(x);
