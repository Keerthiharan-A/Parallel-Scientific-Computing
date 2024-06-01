del = [0.1 0.01 0.005];
red_black = [0.041 16.863 423.849];
serial = [0.005 38.892 780.086];
diag_approach = [0.38 448.157 3628.17];
figure
plot(del, red_black, 'r-*', 'LineWidth', 2, 'DisplayName', 'Red-Black Coloring Approach')
title("Time taken by Solvers vs deta value")
xlabel("Delta")
ylabel("Time taken in s")
hold on
plot(del, serial, 'b--o', 'LineWidth', 2, 'DisplayName', 'Serial Gauss-Seidel')
hold on
plot(del, diag_approach, 'g:s', 'LineWidth', 2, 'DisplayName', 'Diagonal Approach')
legend()
hold off