ngangs = [10, 100, 1000];
time = [0.073282, 0.060716, 0.070454];
figure
plot(ngangs, time, 'b*--', 'LineWidth', 1.5)
xlabel('Number of gangs');
ylabel('Time taken (in s) for complete execution');
title('Time vs No.of gangs plot');