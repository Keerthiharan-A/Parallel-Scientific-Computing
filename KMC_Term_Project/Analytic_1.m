clc
clear
x0 = 0.25;
x_range = 0:0.01:1;
t_values = [0.02, 0.04, 0.08];
T_sums = zeros(length(t_values), length(x_range));
for j = 1:length(t_values)
    t = t_values(j);
    T_sum = zeros(size(x_range));
    for i = 1:numel(x_range)
        x = x_range(i);
        term_sum = 0;
        for n = -5:5
            term = (1/(2*sqrt(pi*t))) * (exp(-(2*n+x-x0)^2/(4*t))-exp(-(2*n+x+x0)^2/(4*t)));
            term_sum = term_sum + term;
        end
        T_sum(i) = term_sum;
    end
    T_sums(j, :) = T_sum;
end
hold on;
for j = 1:length(t_values)
    plot(x_range, T_sums(j, :), 'DisplayName', ['t = ', num2str(t_values(j))]);
end
hold off
xlabel('x');
ylabel('Temperature distribution');
title('T(x,t) vs x');
legend('Location', 'best');