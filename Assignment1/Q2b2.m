thread = [2 4 8];
time_taken = [4.99988 6.99997 3.99995];
figure
plot(thread,time_taken,'b--o',LineWidth=2)
ylim([0,10])
title("Thread vs Time taken plot")
xlabel("Number of threads")
ylabel("Time taken to execute in ms")