clear;clc;
L=10;
EI=3.5e-3;
GIp=3.5e-3;
N=400;
Steps=16000;
Tim=0.2;
Dtim=Tim/Steps;
p=1.0/2;q=1.0/6; r=0;
MN=2.6e6;

% XtimAll=readmatrix('X_Time_All.txt');

Time = readmatrix('Time_All.txt');
Kapa = readmatrix('Kapa_All.txt');
Tau = readmatrix('Tau_All.txt');
Wm=readmatrix('Wm_All.txt');
Wt=readmatrix('Wt_All.txt');
Force_ext=readmatrix('Force_ext.txt');
Real_V=readmatrix('Real_V.txt');

PhaseV = readmatrix('PhaseV_All.txt');
s_int = readmatrix('s_int_All.txt');
Derta = readmatrix('Derta_All.txt');

MKN=3*MN.*Derta;%界面位置宏观迁移率
MKD=48/(pi^2-4)*MN.*Derta;%界面宽度宏观迁移率

Force_ext_1=zeros(Steps,N+1);
Force_ext_2=zeros(Steps,N+1);
Force_ext_3=zeros(Steps,N+1);
Force_ext_Total=zeros(Steps,N+1);

Force_ext_4=zeros(Steps,N);

for i=1:Steps
    for j=1:N+1
        Force_ext_1(i,j)=Force_ext((4*N+3)*(i-1)+4*j-3);
        Force_ext_2(i,j)=Force_ext((4*N+3)*(i-1)+4*j-2);
        Force_ext_3(i,j)=Force_ext((4*N+3)*(i-1)+4*j-1);
        % if j<=N
        %     Fds=XtimAll(407*(i-1)+j+1,3)-XtimAll(407*(i-1)+j,3);
        % end
        % Force=Force_ext_1(i,j)+Force_ext_2(i,j)+Force_ext_3(i,j);
        % if dot(Force,Fds)>0
        %     Force_ext_Total(i,j)=sqrt(Force_ext_1(i,j)^2+Force_ext_2(i,j)^2+Force_ext_3(i,j)^2);
        % else
        %     Force_ext_Total(i,j)=-sqrt(Force_ext_1(i,j)^2+Force_ext_2(i,j)^2+Force_ext_3(i,j)^2);
        % end
        if j<=N
            Force_ext_4(i,j)=Force_ext((4*N+3)*(i-1)+4*j);
        end
    end
end

[li0,kap1i_0,kap2i_0,taui_0,X0]=shape_free(N,p,q,r,L);
[kap1ii_0,kap2ii_0,tauii_0,XX0]=shape_free2(N,p,q,r,L);

kapi_0=sqrt(kap1i_0.^2+kap2i_0.^2);
kapii_0=sqrt(kap1ii_0.^2+kap2ii_0.^2);

n=N-1;
seg=(1:2:2*n-1)/(2*n);
seg=[0 seg 1];
seg=seg*L;
phi=zeros(Steps,N+1);
for i=1:Steps
    phi(i,:)=0.5*(1+tanh((seg-s_int(i))./Derta(i)));
end

dkap=kapi_0-kapii_0;
dtau=taui_0-tauii_0;
Drivingforce=zeros(Steps,N-1);
for i=1:Steps
    for j=1:N-1
        kkk=Kapa((N-1)*(i-1)+j);
        ttt=Tau((N-1)*(i-1)+j);
        Drivingforce(i,j)=dkap(j)*EI*(kkk-kapii_0(j)-phi(i,j+1)*dkap(j))+GIp*dtau(j)*(ttt-tauii_0(j)-phi(i,j+1)*dtau(j));
    end
end

phii=zeros(Steps,N+1);
segi=zeros(Steps,N+1);
ssi=zeros(Steps,N+1);
for i=1:Steps
    k=1;
    for j=1:N+1
        if (s_int(i)-3*Derta(i)<=seg(j)) && (s_int(i)+3*Derta(i)>=seg(j))
            segi(i,k)=seg(j)-s_int(i)+3*Derta(i);
            ssi(i,k)=seg(j);
            k=k+1;
        end
    end
    phii(i,:)=0.5.*(1+tanh((ssi(i,:)-s_int(i))./Derta(i)));
end

PV=PhaseV(1:Steps,1);
PPV=PhaseV(1:Steps,2);
Fdiss_m=Wm(1:Steps,1);
W_T_m=Wm(1:Steps,2);
W_M_m=Wm(1:Steps,3);
W_Double_m=Wm(1:Steps,4);
DFdiss_mm=Wm(1:Steps,5);
DFdiss_mt=Wm(1:Steps,6);
Fdiss_t=Wt(1:Steps,1);
W_T_t=Wt(1:Steps,2);
W_M_t=Wt(1:Steps,3);
W_Double_t=Wt(1:Steps,4);
DFdiss_tt=Wt(1:Steps,5);
DFdiss_tm=Wt(1:Steps,6);

PUS=s_int(1:Steps);
Time=Time(1:Steps);
det=zeros(1,Steps);
lamd=zeros(Steps,2);
DFdissmdt=(W_T_m(2:Steps)-W_T_m(1:Steps-1)+W_M_m(2:Steps)-W_M_m(1:Steps-1))./(0.2/Steps);
DFdisstdt=(W_T_t(2:Steps)-W_T_t(1:Steps-1)+W_M_t(2:Steps)-W_M_t(1:Steps-1))./(0.2/Steps);
%%
for i=1:Steps
    Kss=1-Dtim*MKN(i)*DFdiss_mm(i);
    Ksd=-Dtim*MKN(i)/Derta(i)*Fdiss_m(i)-Dtim*MKN(i)*DFdiss_mt(i);
    Kdd=1-Dtim*MKD(i)/Derta(i)*Fdiss_t(i)-Dtim*MKD(i)*DFdiss_tt(i);
    Kds=-Dtim*MKD(i)*DFdiss_tm(i);
    J=[Kss Ksd;
        Kds Kdd];
    det(i)=Kss*Kdd-Ksd*Kds;
    lamd(i,:)=eig(J);
end
%%
plot(Fdiss_t(800:8425),PhaseV(800:8425,1));
hold on
plot(Time(1:8425),-lamd(1:8425,2));
ylim([0 80])
xlim([0 0.12])
%%
%1
[XX,YY]=meshgrid(seg,Time);
surf(XX,YY,phi)
xlim([0 L])
ylim([0 0.12])
xlabel('Arc length $s$ ($\mu$m)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
ylabel('Time $t$ (s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
xticks([0 2 4 6 8 10])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])
shading interp
colorbar
cb = colorbar;
cb.Label.String = 'Degree of transformation \eta';
cb.Label.FontSize = 16;
clim([0 1]);
colormap(jet);
view([90,-90]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 16);
%%
%2
figure;
hold on;
plot(Time(1:9600),PUS(1:9600),'b','LineWidth',2)

P=polyfit(Time(2800:6800),PUS(2800:6800),1);
P_fit=polyval(P,Time(2800:6800));

% plot(Time(2800:6800),P_fit,'k--','LineWidth',2)
axis([0 0.12 0 10]);
%legend({'Date Point','Fitting Curve'}, 'Fontsize', 14, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$s_{\mathrm{int}}$ ($\mathrm{\mu}$m)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
%
figure;
hold on
cmap = jet(256);
for i=1:8426
    for ii=2:N+1
        if phii(i,ii-1)>=0.9&&phii(i,ii)<=0.1
            kk=ii-1;
            break;
        end
    end
    x_curr = segi(i,1:kk);
    z_curr = phii(i,1:kk);
    t_curr = Time(i)*ones(1,kk);
    for j=1:kk-1
        color_idx = round(z_curr(j) * 255) + 1;
        color_idx = min(max(color_idx, 1), 256);
        
        plot3([x_curr(j), x_curr(j+1)], ...
              [t_curr(j), t_curr(j+1)], ...
              [z_curr(j), z_curr(j+1)], ...
              'Color', cmap(color_idx, :), 'LineWidth', 2);
    end
end
xlabel('Interface region $s-s_{\mathrm{int}}+3D$ ($\mu$m)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
ylabel('Time $t$ (s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
xlim([0 1.8])
ylim([0 0.12])
xticks([0 0.3 0.6 0.9 1.2 1.5 1.8])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])
colormap(jet);
colorbar
cb = colorbar;
cb.Label.String = 'Degree of transformation \eta';
cb.Label.FontSize = 16;
clim([0 1]);
colormap(jet);
view([90,-90]);
set(gca, 'Color', [0 0 0.5]);
hold off
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 16);
%%
%3
% figure;
hold on;
plot(Time(1:8420),Derta(1:8420),'b','LineWidth',2)

% P=polyfit(Time(1600:6000),Derta(1600:6000),0);
% P_fit=polyval(P,Time(1600:6000));
% 
% plot(Time(1600:6000),P_fit,'k--','LineWidth',2)
axis([0 0.12 0 0.3]);
% legend({'Date Point','Fitting Curve'}, 'Fontsize', 14, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 14);
ylabel('$D$ ($\mu$m)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 14);
grid on;
box on;
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 14);
%%
%3
[XX,YY]=meshgrid(seg,Time);
surf(XX,YY,Force_ext_3)
xlim([0 L])
ylim([0 0.12])
xlabel('Arc length $s$ ($\mu$m)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
ylabel('Time $t$ (s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
xticks([0 2 4 6 8 10])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])
shading interp
colorbar
cb = colorbar;
cb.Label.String = 'Fluid viscosity force \it \bf f_{v}';
cb.Label.FontSize = 16;
clim([-2e-6 2e-6]);
colormap(jet);
view([90,-90]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 16);












%%
[XX,YY]=meshgrid(seg(2:end-1),Time);
surf(XX,YY,Drivingforce)
xlim([0 L])
ylim([0 0.12])
xlabel('Arc length $s$ ($\mu$m)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
ylabel('Time $t$ (s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
xticks([0 2 4 6 8 10])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])
shading interp
colorbar
cb = colorbar;
cb.Ticks = [-0.03, -0.02, -0.01, 0, 0.01];
cb.TickLabels = {'-3', '-2', '-1', '0', '1'};
dim = cb.Position;
annotation('textbox', ...
    [dim(1), dim(2)+dim(4)+0.01, 0.05, 0.05], ...
    'String', '×10^{-2}', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 16, ...
    'FontName', 'Times New Roman');
cb.Label.String = 'Driving force $-f_{\mathrm{twist}}-f_{\mathrm{bend}}$';
cb.Label.FontSize = 16;
cb.Label.FontName = 'Times New Roman';
cb.Label.Interpreter = 'Latex';
clim([-0.03 0.01]);
colormap(jet);
view([90,-90]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 16);










%%
figure
yyaxis left;
plot(Time(6000:8000),DFdissmdt(6000:8000).*1000,'r','LineWidth',2);
ylabel('$\dot{F_{\mathrm{e}}^{s}}$ (pN/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-60, 60]);

yyaxis right;
plot(Time(6000:8000),DFdisstdt(6000:8000).*1000,'b','LineWidth',2);
ylabel('$\dot{F_{\mathrm{e}}^{D}}$ (pN/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-6, 6]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$\dot{F_{\mathrm{e}}^{s}}$','$\dot{F_{\mathrm{e}}^{D}}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

xlim([0.075,0.1])
xticks([0.08 0.085 0.09 0.095 0.1])
yyaxis left; 
yticks([-60 -30 0 30 60])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
yticks([-6 -3 0 3 6])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);

%%
%4
figure
yyaxis left;
plot(Time(10:8420),PhaseV(10:8420,1),'r','LineWidth',2);
ylabel('$V_s$ ($\mu$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-160, 320]);

yyaxis right;
plot(Time(10:8420),Fdiss_m(10:8420).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{diss}}^{s}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.4, 0.8]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$V_s$','$F_{\mathrm{diss}}^{s}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
yticks([-160 -80 0 80 160 240 320])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);
%%
%4
figure

yyaxis left;
h_v = plot(Time(60:20:8420),PhaseV(60:20:8420,2),'r','LineWidth',2);
ylabel('$V_D$ ($\mu$m/s)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-40,40]);

yyaxis right;
h_F = plot(Time(60:20:8420),Fdiss_t(60:20:8420).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{diss}}^{D}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.02, 0.02]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$V_D$','$F_{\mathrm{diss}}^{D}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
yticks([-40 -20 0 20 40])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
yticks([-0.02 -0.01 0 0.01 0.02])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);








%%
%6
figure

yyaxis left;
h_v = plot(Time(60:8420),W_T_m(60:8420).*1000+W_M_m(60:8420).*1000,'r','LineWidth',2);
ylabel('$F_{\mathrm{e}}^{s}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.2,0.6]);

yyaxis right;
h_t = plot(Time(60:20:8400),W_T_t(60:20:8400).*1000+W_M_t(60:20:8400).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{e}}^{D}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.01, 0.03]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$F_{\mathrm{e}}^{s}$','$F_{\mathrm{e}}^{D}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
yticks([-0.4 -0.2 0 0.2 0.4 0.6])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
yticks([-0.02 -0.01 0 0.01 0.02 0.03])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);
%%
%7
figure

yyaxis left;
h_v = plot(Time(60:8420),W_Double_m(60:8420).*1000,'r','LineWidth',2);
ylabel('$F_{\mathrm{i}}^{s}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.002,0.003]);

yyaxis right;
h_F = plot(Time(60:20:8400),W_Double_t(60:20:8400).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{i}}^{D}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.004, 0.006]);
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$F_{\mathrm{i}}^{s}$','$F_{\mathrm{i}}^{D}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
yticks([-0.002 -0.001 0 0.001 0.002 0.003])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
yticks([-0.004 -0.002 0 0.002 0.004 0.006])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);










%%
%8
figure

yyaxis left;
h_v = plot(Time(60:8420),W_T_m(60:8420).*1000,'r','LineWidth',2);
ylabel('$F_{\mathrm{twist}}^{s}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-1,1]);

yyaxis right;
h_F = plot(Time(60:8420),W_M_m(60:8420).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{bend}}^{s}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-1, 1]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$F_{\mathrm{twist}}^{s}$','$F_{\mathrm{bend}}^{s}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
% yticks([-0.002 -0.001 0 0.001 0.002 0.003])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
% yticks([-0.004 -0.002 0 0.002 0.004 0.006])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);
%%
%9
figure

yyaxis left;
h_v = plot(Time(60:8410),W_T_t(60:8410).*1000,'r','LineWidth',2);
ylabel('$F_{\mathrm{twist}}^{D}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.04,0.04]);

yyaxis right;
h_F = plot(Time(60:8410),W_M_t(60:8410).*1000,'b','LineWidth',2);
ylabel('$F_{\mathrm{bend}}^{D}$ (pN)', 'interpreter', 'latex', 'FontName', 'Times New Roman', 'Fontsize', 16);
ylim([-0.04, 0.04]);

yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

xlabel('$t$ (s)',  'interpreter', 'latex','FontName', 'Times New Roman', 'Fontsize', 16);
legend({'$F_{\mathrm{twist}}^{D}$','$F_{\mathrm{bend}}^{D}$'}, 'Fontsize', 16, 'Interpreter', 'latex', 'EdgeColor', 'w', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

yyaxis left; 
yticks([-0.04 -0.02 0 0.02 0.04])
set(gca, 'YColor', 'r','FontName', 'Times New Roman', 'Fontsize', 16);
yyaxis right; 
yticks([-0.04 -0.02 0 0.02 0.04])
set(gca, 'YColor', 'b','FontName', 'Times New Roman', 'Fontsize', 16);

