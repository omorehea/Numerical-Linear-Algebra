%{
File: least_sq_poly_fit_Ch.m
Author: Owen Morehead
Purpose: Polynomial fitting solution of the least squares problem using:
         Cholesky Decomposition
%}

%Written by Owen Morehead, Feb 20, 2022

%load x,y data points from data file
dd = load('least_squares_data.dat','-ascii');

%load 3rd and 5th degree best fit polynomial coefficient from data files
d3 = load('fittedCurve_3degree_Ch.dat','ascii');
d5 = load('fittedCurve_5degree_Ch.dat','ascii');


figure(1); grid on; hold on;
plot(dd(:,1),dd(:,2),'oblack','MarkerFaceColor', 'black')

x = 0:0.01:1;

%3rd degree polynomial fit:
%N = 3;
f3 = d3(1) + d3(2)*x + d3(3)*x.^2 + d3(4)*x.^3;

%5th degree polynomial fit
%N = 5;
f5 = d5(1) + d5(2)*x + d5(3)*x.^2 + d5(4)*x.^3 + d5(5)*x.^4 + d5(6)*x.^5;

%exploring 10 and 15 degree polynomial fits. Fortran code exports .dat
%files s.t. the 3rd and 5th degree data file gets overwritten with any new
%polynomial fit. Therefore use d3 and d5 data files.

%f10 = d3(1) + d3(2)*x + d3(3)*x.^2 + d3(4)*x.^3 + d3(5)*x.^4 + ... 
%d3(6)*x.^5 + d3(7)*x.^6 + d3(8)*x.^7 + d3(9)*x.^8 + d3(10)*x.^9 + d3(11)*x.^10;

%f15 = d5(1) + d5(2)*x + d5(3)*x.^2 + d5(4)*x.^3 + d5(5)*x.^4 + ... 
%d5(6)*x.^5 + d5(7)*x.^6 + d5(8)*x.^7 + d5(9)*x.^8 + d5(10)*x.^9 + d5(11)*x.^10 + ...
%d5(12)*x.^11 + d5(13)*x.^12 + d5(14)*x.^13 + d5(15)*x.^14 + d5(16)*x.^15;


plot(x,f3,'r')
plot(x,f5,'b')

xlabel('x', 'Fontsize', 16, 'Interpreter', 'latex')
ylabel('y', 'Fontsize', 16, 'Interpreter', 'latex')
title(sprintf('Cholesky Decomposition Method: Polynomial Fitting Least Squares Problem'), 'Fontsize', 16, 'Interpreter', 'latex')
legend('data','3rd degree polynomial','5th degree polynomial');
