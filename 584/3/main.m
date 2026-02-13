%main part1

close all
clear all
clc


%% Part a
%computation

x_0 = zeros(3,1); % [phi; theta; psi]
t = 0:0.01:10;
opt = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-6);
u = @(t) omega_D_A_D(t);    %omega is expressed as the input

[t_a, x_a] = ode45(@(t,x)funcPartA(t,x,u), t, x_0, opt);
x_a = x_a';
t_a = t_a';
%compute O_mat and evolution of omega as a vector
O_mat_a = O_D_A(x_a(1,:),x_a(2,:),x_a(3,:));
omegas = omega_D_A_D(t_a);

% Plots
figure;
plot(t_a, x_a, 'Linewidth', 2);
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('Euler angles $\vec\alpha$ (rad)','interpreter','latex')
legend({'$\phi$','$\theta$','$\psi$'},'interpreter','latex');
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_a1',[1 0 0]);

figure;
plot(t_a, omegas(1,:), 'Linewidth', 2.5);hold on
plot(t_a, omegas(2,:),'r-.', 'Linewidth', 1.5);hold on
plot(t_a, omegas(3,:), 'Linewidth', 2);hold on
grid on
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('$\omega(t)$ (rad/s)','interpreter','latex')
legend({'$\omega_1$','$\omega_2$','$\omega_3$'},'interpreter','latex');
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_a2',[1 0 0]);


figure;
for ir = 1:9
    r = max(ceil(ir/3),1);
    c = ir- (r-1)*3;
    subplot(3,3,ir);
    plot(t_a, O_mat_a(r,c:3:end), 'Linewidth', 2)
    grid on
    set(gca,'FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ (s)','interpreter','latex')
    ylabel(sprintf('$O_{%d,%d}$',r,c),'interpreter','latex')
end
pos_base = get(gcf,'Position');
set(gcf,'Position',pos_base.*[1 1 2 2]);
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_a3',[1 0 0]);
%% part b
%computation
x_0 = eye(3);
x_0 = reshape(x_0,9,1);
t = 0:0.01:10;
opt = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-6);

u = @(t) omega_D_A_D(t);
%x_b is directly the O_mat evolution expressed as a 9x1*length(t) vector
[t_b, x_b] = ode45(@(t,x)poisson(t,x,u), t, x_0, opt);
x_b = x_b';
t_b = t_b';
%compute omegas
omegas = omega_D_A_D(t_a);
%reshape the O_mat to get 3x3*length(t)
O_mat_b = reshape(x_b,3,[]);

%plot
figure;
for ir = 1:9
    r = max(ceil(ir/3),1);
    c = ir- (r-1)*3;
    subplot(3,3,ir);
    plot(t_a, O_mat_b(r,c:3:end), 'Linewidth', 2);hold on;
    plot(t_a, O_mat_a(r,c:3:end),'-.', 'Linewidth', 1.5);
    
    grid on
    set(gca,'FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$t$ (s)','interpreter','latex')
    ylabel(sprintf('$O_{%d,%d}$',r,c),'interpreter','latex')
%     legend({'$\phi$','$\theta$','$\psi$'},'interpreter','latex');
end
legend({'Poisson $O_{D/A}$','Euler $O_{D/A}$'}, 'interpreter','latex');

pos_base = get(gcf,'Position');
set(gcf,'Position',pos_base.*[1 1 2 2]);
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_b1',[1 0 0]);


%% part C
%computation

ang_b = angleFrom_O_mat(O_mat_b);
ang_a = angleFrom_O_mat(O_mat_a);


%plot
figure;
h1 = plot(t_a, x_a,'b-', 'Linewidth', 3);hold on;
h2 = plot(t_a, ang_a,'b--','Color','#D95319', 'Linewidth', 2.5);hold on;
h3 = plot(t_b, ang_b,'y-.', 'Linewidth', 1);hold on;


grid on
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('angle (rad)','interpreter','latex')
legend([h1(1),h2(1),h3(1)],{'from a)','from $O_{D/A}$ of a)','from $O_{D/A}$ of b)'},'interpreter','latex');
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_c1',[1 0 0]);

figure;
subplot(3,1,1);
plot(t_a,x_a-ang_a);
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('$\Delta(\vec\alpha)$','interpreter','latex')
title('$\Delta$ E.A. = E.A. from a) - E.A. from $O_{D/A,a}$','interpreter','latex')

subplot(3,1,2);
plot(t_b,x_a-ang_b)
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('$\Delta(\vec\alpha)$','interpreter','latex')
title('$\Delta$ E.A. = E.A. from a) - E.A. from $O_{D/A,b}$','interpreter','latex')

subplot(3,1,3);
plot(t_a,ang_a-ang_b);
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$ (s)','interpreter','latex')
ylabel('$\Delta(\vec\alpha)$','interpreter','latex')
title('  $\Delta$ E.A. = E.A. from $O_{D/A,a}$ - E.A. from $O_{D/A,b}$','interpreter','latex')


pos_base = get(gcf,'Position');
set(gcf,'Position',pos_base.*[1 1 1 1.4]);
%save figures
fg_saveAs_pdf_png_fig([],'Problem7_c2',[1 0 0]);



