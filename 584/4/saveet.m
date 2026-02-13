save('euler_01.mat','O')

dt = 0.1;
timeVec = [0:0.1:10];
0 = squeese(O);

figure(1)
for ir = 1:9
    r = max(ceil(ir/3),1);
    c = ir- (r-1)*3;
    subplot(3,3,ir);
    plot(timeVec, O(r,c:3:end),'LineWidth',1.5,'Color',[0, 0.5, 0.5]);
    grid on
    xlabel('time [s]')
    ylabel(sprintf('O_{%d,%d}',r,c))
end
legend('O_{B/E}(t)','Location','best');
sgtitle('Problem 3: Components of O_{B/E} vs Time')