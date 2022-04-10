%{ 
File: iteration_algs_errors.m
Author: Owen Morehead
Purpose: Reads and plots all error data as .dat files outputted
form Fortran program. Error is calculated at each iteration for
Gauss-Jacobi and Gauss-Seidel iterative algorithms.
%}


%------ Plots of error for Gauss-Jacobi and Gauss-Seidel Methods -----
%------ data output from Fortran code -----
%both GJ and GS converge for: D = 10,100,1000

dS10 = load('errorS_100010.dat','-ascii');
dS100 = load('errorS_100100.dat','-ascii');
dS1000 = load('errorS_101000.dat','-ascii');

dJ10 = load('errorJ_100010.dat','-ascii');
dJ100 = load('errorJ_100100.dat','-ascii');
dJ1000 = load('errorJ_101000.dat','-ascii');



figure(1); grid on; hold on;
set(gca, 'yscale', 'log')

plot(1:length(dS10),dS10(:),'o-b','MarkerFaceColor', 'b')
plot(1:length(dJ10),dJ10(:),'o-r','MarkerFaceColor', 'r')

xlabel('Iteration', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('log[Error: $||b - Ax||_2]$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Convergence Error for Iterative Methods','D = 10'}, 'Fontsize', 20, 'Interpreter', 'latex')
legend('Gauss-Seidel','Gauss-Jacobi','fontsize',18,'interpreter','latex');

figure(2); grid on; hold on;
set(gca, 'yscale', 'log')

plot(1:length(dS100),dS100(:),'o-b','MarkerFaceColor', 'b')
plot(1:length(dJ100),dJ100(:),'o-r','MarkerFaceColor', 'r')

xlabel('Iteration', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('log[Error: $||b - Ax||_2]$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Convergence Error for Iterative Methods','D = 100'}, 'Fontsize', 20, 'Interpreter', 'latex')
legend('Gauss-Seidel','Gauss-Jacobi','fontsize',18,'interpreter','latex');


figure(3); grid on; hold on;
set(gca, 'yscale', 'log')

plot(1:length(dS1000),dS1000(:),'o-b','MarkerFaceColor', 'b')
plot(1:length(dJ1000),dJ1000(:),'o-r','MarkerFaceColor', 'r')

xlabel('Iteration', 'Fontsize', 20, 'Interpreter', 'latex')
ylabel('log[Error: $||b - Ax||_2]$', 'Fontsize', 20, 'Interpreter', 'latex')
title({'Convergence Error for Iterative Methods','D = 1000'}, 'Fontsize', 20, 'Interpreter', 'latex')
legend('Gauss-Seidel','Gauss-Jacobi','fontsize',18,'interpreter','latex');

