%AM213A Problem 12
%Written by Owen Morehead, Jan 23 2021

%Comparing function with two definitions
%function ill conditioned at x = 2

x = 1.92:0.001:2.080;
f = (x-2).^9;
g = x.^9 - 18*x.^8 + 144*x.^7 - 672*x.^6 + 2016*x.^5 - 4032*x.^4 ...
+ 5376*x.^3 - 4608*x.^2 + 2304*x - 512;

figure(1); grid on; hold on;
scatter(x,f,'r');
scatter(x,g,'b');

xlabel(sprintf("x"),'Interpreter','latex', 'fontsize',14);
%ylabel(sprintf("Functions F(x) and G(x)"),'Interpreter','latex','fontsize',14);
legend({sprintf('$f(x) = (x-2)^9$'),sprintf('$g(x) = x^9 - 18x^8 + 144x^7 - 672x^6 + 2016x^5 - 4032x^4 + 5376x^3 - 4608x^2 + 2304x - 512$')},'fontsize',11,'Interpreter','latex')
title("Problem 12: f(x) vs. g(x)",'Interpreter','latex','fontsize',16);

