%AM213A HW2 Problem 5 
%Equation of plane found through three points
%Written by Owen Morehead, Feb 9 2022

clc
clear all
close all

%different formats for the three points (x,y,z) 
p1 = [1,2,3]; p2 = [-3,2,5]; p3  = [pi,exp(1),-sqrt(2)];
pxs = [1,-3,pi]; pys = [2,2,exp(1)]; pzs = [3,5,-sqrt(2)];

[x, y] = meshgrid(-8:.1:8);

%%below first 2 equivalent plane eq, z, for a,b,c not normalized

%z = ((3.9034e-2*(x-pxs(3)) + 0.3634*(y-pys(3)))/-7.8067e-2) + pzs(3);
%z = (3.9034e-2*x + 0.3634*y - 1)/(-7.8067e-2); %equivalent definition

%equivalent definition different d = ..  -> normalized wrt a
%can use coordinates of any point (x_0,y_0,z_0) in plane equaiton

%z = (1*(x) + 9.3095*(y) - 25.6189)/(-2.0); 
z = ((1*(x-pxs(3)) + 9.3095*(y-pys(3)))/-2) + pzs(3);
surf(x,y,z,'FaceAlpha',0.5,'EdgeColor','none','FaceColor',[.1,.5,.9]);

hold on;

scatter3(p1(1),p1(2),p1(3),1000,'r.');
scatter3(p2(1),p2(2),p2(3),1000,'b.');
scatter3(p3(1),p3(2),p3(3),1000,'g.');

xlabel('x','Interpreter','latex', 'fontsize',14); 
ylabel('y','Interpreter','latex', 'fontsize',14); 
zlabel('z','Interpreter','latex', 'fontsize',14)
title("Plane Equation: $z = -\frac{a(x-x_0) + b(y-y_0)}{c} + z_0$",'Interpreter','latex', 'fontsize',16) 
legend({sprintf(''),sprintf('A = (%.5f,%.5f,%.5f)',p1(1),p1(2),p1(3)),sprintf('B = (%.5f,%.5f,%.5f)',p2(1),p2(2),p2(3)),sprintf('C = (%.5f,%.5f,%.5f)',p3(1),p3(2),p3(3))},'fontsize',15,'Interpreter','latex')
