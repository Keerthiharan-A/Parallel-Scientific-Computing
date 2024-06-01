clc
clear
y1 = [5.10572 4.09424 1.81829 -1.13726 -3.68363 -4.94639 -4.48038 -2.44947 0.437169 3.17108 4.79724 4.74759 3.03947 0.269574 -2.59449 -4.55223 -4.91974 -3.56865 -0.970901 1.96591 4.21636 4.99245 4.02999 1.63939 -1.24794 -3.98281];
x1 = linspace(0,3,26);
x2 = linspace(0,3,100);
y2 = 5*cos(5*x2);
figure
scatter(x1, y1, 'r*', 'LineWidth', 2, 'DisplayName', 'LU Decomposition solution')
xlabel("X")
ylabel("Y")
hold on
plot(x2, y2, 'b-', 'LineWidth', 1, 'DisplayName', 'Exact Solution')
legend()