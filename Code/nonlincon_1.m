function [c ceq] = nonlincon_1(t)
syms x y %  L r rp h b
% r = Radius of /moving platform
% rp = Radius of hole in moving platform
% h = Height(Pad) of link
% b = Width of link

% L = 0.02;
% r = 0.015;
% rp = 0.002;
% h = 0.010;
% b = 0.010;

B1x = -0.0125;
B1y = -0.0216505;
P1x = -0.009802;
P1y = 0.001978;

Eq_1_1 = (x - B1x)^2 + (y - B1y)^2 == t(1)^2;
Eq_1_2 = (x - P1x)^2 + (y - P1y)^2 == t(1)^2;

Solution_1 = solve([Eq_1_1,Eq_1_2],[x,y]);
xSolution_1 = vpa(Solution_1.x,5) %solve(Eq_1_1,x);%L1x
ySolution_1 = vpa(Solution_1.y,5) %solve(Eq_1_2,y);%L1y

% Alpha_1 =atan((ySolution_1+B1y)/(xSolution_1-B1x));
% Alpha_1_Degree = rad2deg(Alpha_1)

Gamma_1 = acos((B1x-xSolution_1)/t(1));
Gamma_1_Degree = rad2deg(Gamma_1);
Alpha_1_Degree = vpa(180-Gamma_1_Degree,5)

B1P1=sqrt((B1x-P1x)^2+(B1y-P1y)^2);

Theta_1 =acos((2*t(1)^2-B1P1^2)/2*t(1)^2);
Theta_1_Degree = rad2deg(Theta_1);

Beta_1 = vpa(180+Theta_1_Degree,5)


B2x = 0;
B2y = 0.025;
P2x = 0.006614;
P2y = 0.00750;
Eq_2_1 = (x - B2x)^2 + (y - B2y)^2 == t(1)^2;
Eq_2_2 = (x - P2x)^2 + (y - P2y)^2 == t(1)^2;
 
Solution_2 = solve([Eq_2_1,Eq_2_2],[x,y]);
xSolution_2 = vpa(Solution_2.x,5) % solve(Eq_2_1,x);
ySolution_2 = vpa(Solution_2.y,5) % solve(Eq_2_2,y);

% Alpha_2 =atan((ySolution_2+B2y)/(xSolution_2-B2x));
% Alpha_2_Degree = rad2deg(Alpha_2)

Gamma_2 = asin((xSolution_2-B2x)/t(1));
Gamma_2_Degree = rad2deg(Gamma_2);
Alpha_2_Degree = vpa(270+Gamma_2_Degree,5)

B2P2=sqrt((B2x-P2x)^2+(B2y-P2y)^2);

Theta_2 =acos((2*t(1)^2-B2P2^2)/2*t(1)^2);
Theta_2_Degree = rad2deg(Theta_2);

Beta_2 = vpa(180+Theta_2_Degree,5)


B3x = 0.0125;
B3y = -0.0216505;
P3x = 0.003188;
P3y = -0.009478;
Eq_3_1 = (x - B3x)^2 + (y - B3y)^2 == t(1)^2;4
Eq_3_2 = (x - P3x)^2 + (y - P3y)^2 == t(1)^2;

Solution_3 = solve([Eq_3_1,Eq_3_2],[x,y]);
xSolution_3 = vpa(Solution_3.x,5) % solve(Eq_3_1,x);
ySolution_3 = vpa(Solution_3.y,5) % solve(Eq_3_2,y);
 
% Alpha_3 =atan((ySolution_3+B3y)/(xSolution_3-B3x));
% Alpha_3_Degree = rad2deg(Alpha_3)

Gamma_3 = acos((B3x-xSolution_3)/t(1));
Gamma_3_Degree = rad2deg(Gamma_3);
Alpha_3_Degree = vpa(Gamma_3_Degree,5)

B3P3=sqrt((B3x-P3x)^2+(B3y-P3y)^2);

Theta_3 =acos((2*t(1)^2-B3P3^2)/2*t(1)^2);
Theta_3_Degree = rad2deg(Theta_3);

Beta_3 = vpa(180-Theta_3_Degree,5)

 
Forward_J = [cos(Alpha_1_Degree(1)+Beta_1) sin(Alpha_1_Degree(1)+Beta_1) (t(2)/2);cos(Alpha_2_Degree(1)+Beta_2) sin(Alpha_2_Degree(1)+Beta_2) (t(2)/2);cos(Alpha_3_Degree(1)+Beta_3) sin(Alpha_3_Degree(1)+Beta_3) (t(2)/2)]
Inverse_J = [cos(Alpha_1_Degree(2)+Beta_1) sin(Alpha_1_Degree(2)+Beta_1) (t(2)/2);cos(Alpha_2_Degree(2)+Beta_2) sin(Alpha_2_Degree(2)+Beta_2) (t(2)/2);cos(Alpha_3_Degree(2)+Beta_3) sin(Alpha_3_Degree(2)+Beta_3) (t(2)/2)]
 
A_Inverse = inv(Forward_J)
B_Inverse = inv(Inverse_J)
 
E = [0 -1;1 0]
B01 = t(1)*[sin(Beta_1) 0 0;0 sin(Beta_2) 0;0 0 sin(Beta_3)]
 
J_1 = A_Inverse*B01
J_2 = B_Inverse*B01
 
J_Inverse_1 = inv(J_1)
J_Inverse_2 = inv(J_2)
 
% Size_of_J_Inverse_1=size(J_Inverse_1)
 
% Size_of_J_1=size(J_1)
 
% Size_of_J_Inverse_2=size(J_Inverse_2)
 
% Size_of_J_2=size(J_2)
 
Condition_Num_1 = vpa(1/(norm(J_1)*norm(J_Inverse_1)),5)
Condition_Num_2 = vpa(1/(norm(J_2)*norm(J_Inverse_2)),5)

Shear_Mod= 7.93*10^10;
Youngs_Mod= 210*10^9;
Iz = (t(4)*t(5)^3)/12;
Iy = (t(4)^3)/12;
Ix = Iy+Iz;
ALb = t(4)*t(5);
K_Leg_Inv = [(t(1))/(Youngs_Mod*ALb) 0 0 0 0 0;0 (t(1))/(3*Iz) 0 0 0 (t(1)^2)/(2*Youngs_Mod*Iz);0 0 (t(1)^3)/(3*Youngs_Mod*Iy) 0 (-t(1)^2)/(2*Youngs_Mod*Iy) 0;0 0 0 (t(1))/(Shear_Mod*Ix) 0 0;0 0 (-t(1)^2)/(2*Youngs_Mod*Iy) 0 (t(1))/(Youngs_Mod*Iy) 0;0 (t(1)^2)/(2*Youngs_Mod*Iz) 0 0 0 (t(1))/(Youngs_Mod*Iz)];
%Stifness matrix of intermideate legs
K_MP_Inv = [(t(2))/(Youngs_Mod*ALb) 0 0 0 0 0;0 (t(2))/(3*Iz) 0 0 0 (t(2)^2)/(2*Youngs_Mod*Iz);0 0 (t(2)^3)/(3*Youngs_Mod*Iy) 0 (-t(2)^2)/(2*Youngs_Mod*Iy) 0;0 0 0 (t(2))/(Shear_Mod*Ix) 0 0;0 0 (-t(2)^2)/(2*Youngs_Mod*Iy) 0 (t(2))/(Youngs_Mod*Iy) 0;0 (t(2)^2)/(2*Youngs_Mod*Iz) 0 0 0 (t(2))/(Youngs_Mod*Iz)];
%Stifness matrix of Moving Platform
K_Act_Inv =1/(t(1)/(Youngs_Mod*ALb));
%Stifness matrix of Actuator
K_Theta = blkdiag(K_Act_Inv,K_Leg_Inv,K_MP_Inv)
%Combining all 3 matrices into one diagonal matrix
K_Theta_Inv = inv(K_Theta)


One_Matrix = ones(3,26);
J_1_One_Mat = J_1*One_Matrix;
J_1_One_Mat_Reshape = reshape(J_1_One_Mat,13,6)
%Converting 3x26 matrix into 13x6
Trans_of_J_1_One_Mat_Reshape = transpose(J_1_One_Mat_Reshape);
J_2_One_Mat = J_2*One_Matrix;
J_2_One_Mat_Reshape = reshape(J_2_One_Mat,13,6)
Trans_of_J_2_One_Mat_Reshape = transpose(J_2_One_Mat_Reshape);

Stiffness_1 = 3*Trans_of_J_1_One_Mat_Reshape*K_Theta_Inv*J_1_One_Mat_Reshape
%Norm_of_Stiffness_1 = vpa(norm(Stiffness_1),2)
Stiffness_2 = 3*Trans_of_J_2_One_Mat_Reshape*K_Theta_Inv*J_2_One_Mat_Reshape
%Norm_of_Stiffness_2 = vpa(norm(Stiffness_2),2)

 c(1)= -double(Condition_Num_1)+0.1; %inverse conditioning number1
 c(2)= -double(Condition_Num_2)+0.1; %inverse conditioning number2
 c(3) = double(sqrt(2*t(1))^2-((sqrt(3)/2*(t(5)-t(2)))^2)+(t(2)/2)-(t(5)/2)-50)
  % c(5)= -Norm_of_Stiffness_1+10^5; %condition of force using stiffness matrix1 %stiff_norm1
% c(4)= -Norm_of_Stiffness_2+10^5; %condition of force using stiffness matrix2 %stiff_norm2

d = 7850;%Density
Mass_of_MP = pi*(t(2)^2)*t(4)*d-3*pi*(t(3)^2)*t(4)*d;
Mass_of_Links = t(1)*t(5)*t(4)*d ;
Total_Mass = Mass_of_MP + 6*Mass_of_Links;
ceq=[];

end

