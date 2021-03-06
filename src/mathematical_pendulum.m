f = @(x,y) [y(2);-sin(y(1))];
x0 = 0;
xN = 4 * 1/sqrt(2) * integral(@(theta) 1./sqrt(cos(theta) - cos(pi/3)),0,pi/3);
y0 = [-pi/3;0];
[x,y] = DOPRI5(f,x0,xN,y0);
plot(y(1,:),y(2,:),'-','color','black');
grid on;
l = legend('\texttt{DOPRI5}','location','southeast');
set(l, 'fontsize', 18, 'interpreter', 'latex');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0, 0, 12, 6];
set(fig, 'PaperSize', [12, 6]);
saveas(fig, 'mathematical_pendulum.pdf');