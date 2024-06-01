nt = [2 4 8 16];
red_black = [447.128 415.961 423.849 432.303];
diag_solver = [2335.63 3155.19 3684.42 4084.42];
figure
plot(nt, red_black, 'r-*', 'LineWidth', 2, 'DisplayName', 'Red-Black Coloring Approach')
title("Time taken by Solvers vs number of threads")
xlabel("Thread count")
ylabel("Time taken in s")
hold on
plot(nt, diag_solver, 'g:s', 'LineWidth', 2, 'DisplayName', 'Diagonal Approach')
legend()
hold off