red_balck = [0 0.1414 0.267854 0.379386 0.476018 0.557769 0.624654 0.676684 0.713869 0.736211 0.743711 0.736365 0.714165 0.677101 0.625158 0.55832 0.476568 0.379882 0.26824 0.14162 0];
y2 = [-0 0.1425 0.27 0.3825 0.48 0.5625 0.63 0.6825 0.72 0.7425 0.75 0.7425 0.72 0.6825 0.63 0.5625 0.48 0.3825 0.27 0.1425 -0];
diag_method = [0 0.1414 0.267854 0.379386 0.476018 0.557769 0.624654 0.676684 0.713869 0.736211 0.743711 0.736365 0.714165 0.677101 0.625158 0.55832 0.476568 0.379882 0.26824 0.14162 0];
x = linspace(-1,1,21);
figure
plot(x, red_balck, 'r-', 'LineWidth', 4, 'DisplayName', 'Serial Gauss Seidel Solution')
xlabel("X")
ylabel("Y")
hold on
plot(x, y2, 'b-', 'LineWidth', 2, 'DisplayName', 'Exact Solution')
hold on
plot(x, diag_method, 'g--', 'LineWidth', 3, 'DisplayName', 'Diagonal approach')
legend()