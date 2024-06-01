clc
clear
jacobi_soln_1 = [0 0.0339049 0.059572 0.0787079 0.0929142 0.103648 0.112197 0.119657 0.126932 0.134725 0.14354 0.153686 0.165275 0.178222 0.192245 0.206856 0.22136 0.234839 0.24615 0.253912 0.256499];
jacobi_soln_2 = [0 0.0402911 0.0714942 0.094671 0.111112 0.122304 0.129778 0.1349 0.138662 0.141549 0.14354 0.144223 0.142988 0.139226 0.132451 0.122304 0.108438 0.0903445 0.0671678 0.0376173 0];
% It took 716 iterations to converge to the error 9.93947e-05
x = linspace(-1,1,21);
figure
plot(x, jacobi_soln_1, 'r*--', 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('phi at y = 0');
title('Phi at y = 0 vs x');
legend('Jacobi Solution');
figure
plot(x, jacobi_soln_2, 'b*--', 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('phi at x = 0');
title('Phi at x = 0 vs y');
legend('Jacobi Solution');