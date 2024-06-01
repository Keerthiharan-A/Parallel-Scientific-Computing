serial = [1e-6 0.003 1.157];
parallel = [0.0051750 0.0068163 1.031320];
n = [10 100 1000];
plot(n,serial,'b*--', 'LineWidth', 1, 'DisplayName', 'Serial solver')
hold on
plot(n,parallel,'r*--', 'LineWidth', 1, 'DisplayName', 'Parallel solver')
xlabel('Matrix size')
ylabel('Time taken (in s)')
title('Comparing Time taken in Serial and Parallel execution')
legend()