clear;clc;
t=50:10:150;
p1=[78 94 108 122 135 148 161 173 186 197 209];
p2=[84 101 109 117 135 141 156 160 170 189 194];
p3=[85 93 105 122 133 147 155 170 175 190 202];
figure;
hold on;

P1=polyfit(t,p1,1);
P_fit1=polyval(P1,t);

P2=polyfit(t,p2,1);
P_fit2=polyval(P2,t);

P3=polyfit(t,p3,1);
P_fit3=polyval(P3,t);

scatter(t, p1,80, 'blue');
plot(t,P_fit1,'b--','LineWidth',1)
scatter(t, p2,80, 'black');
plot(t,P_fit2,'k--','LineWidth',1)
scatter(t, p3,80, 'red','filled');
plot(t,P_fit3,'r','LineWidth',1)
axis([40 160 50 250]);
legend({'Simulation data (Ko et al.)','Fitting curve (Ko et al.)',...
    'Simulation data (Xu et al.)','Fitting curve (Xu et al.)',...
    'Simulation data (This work)','Fitting curve (This work)'},...
    'Fontsize', 12, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
xlabel('$f$ (Hz)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mu$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
clear;clc;
x=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
y1=[0 76 104 122 133 136 138 141 143 146 147];
y2=[0 79 103 120 128 132 135 138 141 144 146];

ft = fittype('a*log(b*x+1)');
opts1 = fitoptions(ft);
opts1.StartPoint = [70, 20];
fitted_model1 = fit(x', y1', ft,opts1);

opts2 = fitoptions(ft);
opts2.StartPoint = [70, 20];
fitted_model2 = fit(x', y2', ft,opts2);

x_fit=0:0.01:1;
y1_fit=fitted_model1(x_fit);
y2_fit=fitted_model2(x_fit);

figure;
hold on;
scatter(x, y1,80, 'blue');
plot(x_fit',y1_fit,'b--','LineWidth',1);
scatter(x, y2,80, 'red','filled');
plot(x_fit',y2_fit,'r','LineWidth',1);

axis([0 1.1 0 200]);
legend({'Simulation data (Xu et al.)','Fitting curve (Xu et al.)',...
    'Simulation data (This work)','Fitting curve (This work)'},...
    'Fontsize', 12, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
xlabel('$\mathrm{\lambda}$',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mathrm{\mu}$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);




%%
%100
clear;clc;
x=[2000 2500 3000 3500 4000];
y1=[114 126 136 145 152];
figure;
hold on;
plot(x,y1,'b--','LineWidth',2)

scatter(x, y1,80, 'blue');

axis([1500 4500 100 160]);
xlabel('Time steps',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mathrm{\mu}$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
ax=gca;
ax.XTick=[2000 2500 3000 3500 4000];
ax.XTickLabel={'2','2.5','3','3.5','4'};

xlims=ax.XLim;
ylims=ax.YLim;

text(xlims(2),ylims(1)-(ylims(2)-ylims(1))*0.03,'\times10^3','FontSize',12,'HorizontalAlignment','left');
yticks([80 100 120 140 160])
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
%134
clear;clc;
x=[2000 3000 4000 5000 6000];
y1=[104 131 142 148 155];
figure;
hold on;
plot(x,y1,'b--','LineWidth',2)

scatter(x, y1,80, 'blue');

axis([1500 6500 100 160]);
xlabel('Time steps',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mathrm{\mu}$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
ax=gca;
ax.XTick=[2000 3000 4000 5000 6000];
ax.XTickLabel={'2','3','4','5','6'};

xlims=ax.XLim;
ylims=ax.YLim;

text(xlims(2),ylims(1)-(ylims(2)-ylims(1))*0.03,'\times10^3','FontSize',12,'HorizontalAlignment','left');
yticks([80 100 120 140 160])
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
%200
clear;clc;
x=[2000 4000 6000 8000 10000];
y1=[86 125 148 149 150];
figure;
hold on;
plot(x,y1,'b--','LineWidth',2)

scatter(x, y1,80, 'blue');

axis([1000 11000 80 160]);
%legend({'Phase field model','Two parameter model'}, 'Fontsize', 14, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
xlabel('Time steps',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mathrm{\mu}$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
ax=gca;
ax.XTick=[2000 4000 6000 8000 10000];
ax.XTickLabel={'2','4','6','8','10'};

xlims=ax.XLim;
ylims=ax.YLim;

text(xlims(2),ylims(1)-(ylims(2)-ylims(1))*0.03,'\times10^3','FontSize',12,'HorizontalAlignment','left');
yticks([80 100 120 140 160])
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
%400
clear;clc;
x=[4000 8000 12000 16000 20000];
y1=[86 126 145 149 150];
figure;
hold on;
plot(x,y1,'b--','LineWidth',2)

scatter(x, y1,80, 'blue');

axis([2000 22000 80 160]);
xlabel('Time steps',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$V_{s}$ ($\mathrm{\mu}$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
ax=gca;
ax.XTick=[4000 8000 12000 16000 20000];
ax.XTickLabel={'4','8','12','16','20'};

xlims=ax.XLim;
ylims=ax.YLim;

text(xlims(2),ylims(1)-(ylims(2)-ylims(1))*0.03,'\times10^3','FontSize',12,'HorizontalAlignment','left');
yticks([80 100 120 140 160])
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
