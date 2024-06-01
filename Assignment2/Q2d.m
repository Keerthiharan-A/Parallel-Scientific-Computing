clc
clear
nproc = [2,4,8,16];
j_tseq = 2078.19;
j_time = [1119.43, 1424.88, 1820.9, 2232.21];
j_speedup = j_tseq./j_time;
j_efficiency = j_speedup./nproc;
gs_tseq = 1052.56;
gs_time = [610.614, 694.991, 828.622, 1056.92];
g_speedup = gs_tseq./gs_time;
g_efficiency = g_speedup./nproc;
figure
plot(nproc, j_speedup, 'r*--', 'LineWidth', 1);
xlabel('Number of processors');
ylabel('Speedup');
xlim([0,16]);
ylim([0.4, 2]);
hold on
plot(nproc, g_speedup, 'b*--', 'LineWidth', 1);
xlabel('Number of processors');
ylabel('Speedup');
legend('Jacobi SpeedUp', 'Red Black SpeedUp')
title('Jacobi and Red Black SpeedUp');

figure
plot(nproc, j_efficiency, 'r*--', 'LineWidth', 1);
xlabel('Number of processors');
ylabel('Efficiency');
xlim([0,16]);
hold on
plot(nproc, g_efficiency, 'b*--', 'LineWidth', 1);
xlabel('Number of processors');
ylabel('Efficiency');
legend('Jacobi Efficiency', 'Red Black Efficiency')
title('Jacobi and Red Black Efficiency');