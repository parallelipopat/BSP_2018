% Plots the Lennard-Jones potential curve

r = 0:0.01:5;
inv_r6 = r.^-6;
inv_r12 = inv_r6.^2;
V = 4*(inv_r12-inv_r6);
plot(r(1,:),V(1,:),'b');
xlim([0 3]);
ylim([-1.5 2]);
hold on;
plot(r,zeros(size(r,2)),'k');
set(gca,'xtick',[],'ytick',[0]);
ylabel('V');
xlabel('r','fontweight','bold');
