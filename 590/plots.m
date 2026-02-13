% plot the orbit
if FlagOrbitPlot == 1
    figure(1); % 3D Position
    plot3(Xorb(:,1), Xorb(:,2), Xorb(:,3),'g-','linewidth',3);
    xlabel('x (km)','fontsize',14);ylabel('y (km)','fontsize',14);
    zlabel('z (km)','fontsize',14);
    set(gca,'fontsize',14)
    [XS, YS, ZS] = sphere(30); % plot the Earth
    hold on;
    surf(XS*Re, YS*Re, ZS*Re);
    axis equal
    %     title('3D Trajectory Conventional EOM')
    hold off
    grid on
end
if orbitRep == 1
    plotTimeVec = Torb;
    xLabel = 'Time (s)';
elseif orbitRep > 1
    plotTimeVec = Torb./period;
    xLabel = 'Orbits';
end

if FlagAngVelPlot == 1
    % plot sc angular velocity

    fig = figure(2);
    subplot(3,1,1)
    plot(plotTimeVec,Xorb(:,7),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \omega_x (rad/s)$','interpreter','latex'); grid on;

    subplot(3,1,2)
    plot(plotTimeVec,Xorb(:,8),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \omega_y (rad/s)$','interpreter','latex'); grid on;

    subplot(3,1,3)
    plot(plotTimeVec,Xorb(:,9),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \omega_z (rad/s)$','interpreter','latex'); grid on;
    if FlagControl == 0
        titleName = 'SC Dynamics Angular Velocity No Control';
    elseif FlagControl == 1
        titleName = 'SC Dynamics Angular Velocity Saturated Control';
    elseif FlagControl == 2
        titleName = 'SC Dynamics Angular Velocity Unsaturated Control';
    end
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')
end



if FlagMagFieldPlot == 1 || flagPlotDipMoment == 1
    for ii=1:length(timeVec)
        [~,magField,dipMom] = scDyn13(timeVec(ii),Xorb(ii,:)',TimeStamp0,FlagMagField,FlagControl,constRot);
        if FlagMagFieldPlot == 1
            XYZ_NED = magField(:,1);
            XYZ_ECEF = magField(:,2);
            XYZ_ECI = magField(:,3);
            XYZ_P = magField(:,4);
            magFieldMatrix_NED(ii,:) = XYZ_NED';
            magFieldMatrix_ECEF(ii,:) = XYZ_ECEF';
            magFieldMatrix_ECI(ii,:) = XYZ_ECI';
            magFieldMatrix_P(ii,:) = XYZ_P';
        end
        if flagPlotDipMoment == 1
            dipMomMatrix(ii,:) = dipMom;
        end
    end
end
7
% Magnetic Field Components
if FlagMagFieldPlot == 1
    fig = figure(3); % Magnetic Field Components and Magnitude in NED
    plot(plotTimeVec, magFieldMatrix_NED(:,1), 'Color',[0 .5 .5],'linewidth',3);
    hold on
    plot(plotTimeVec, magFieldMatrix_NED(:,2), 'Color',[.5 0 .5],'linewidth',3);
    plot(plotTimeVec, magFieldMatrix_NED(:,3), 'Color',[0.7 0.6 0.3],'linewidth',3);
    plot(plotTimeVec, vecnorm(magFieldMatrix_NED,2,2), 'Color',[0.1 0.1 0.1],'linewidth',3);
    legend('$ \stackrel{\rightharpoonup}{b}|_{NED,x} $', '$ \stackrel{\rightharpoonup}{b}|_{NED,y} $', '$ \stackrel{\rightharpoonup}{b}|_{NED,z} $', '$ |\stackrel{\rightharpoonup}{b}| $','interpreter','latex','Location','best','FontSize',16)
    titleName = 'Magnetic Field on Orbit Resolved in NED';
    %     title(titleName)
    ylabel('Magnetic Flux Density (nT)','interpreter','latex'); xlabel(xLabel,'interpreter','latex'); set(gca,'fontsize',12); grid on
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')


    fig = figure(4); % Magnetic Field Components and Magnitude in ECEF
    plot(plotTimeVec, magFieldMatrix_ECEF(:,1), 'Color',[0 .5 .5],'linewidth',3);
    hold on
    plot(plotTimeVec, magFieldMatrix_ECEF(:,2), 'Color',[.5 0 .5],'linewidth',3);
    plot(plotTimeVec, magFieldMatrix_ECEF(:,3), 'Color',[0.7 0.6 0.3],'linewidth',3);
    plot(plotTimeVec, vecnorm(magFieldMatrix_ECEF,2,2), 'Color',[0.1 0.1 0.1],'linewidth',3);
    legend('$ \stackrel{\rightharpoonup}{b}|_{E,x} $', '$ \stackrel{\rightharpoonup}{b}|_{E,y} $', '$ \stackrel{\rightharpoonup}{b}|_{E,z} $', '$ |\stackrel{\rightharpoonup}{b}| $','interpreter','latex','Location','best','FontSize',16)
    titleName = 'Magnetic Field on Orbit Resolved in ECEF';
    %     title(titleName)
    ylabel('Magnetic Flux Density (nT)','interpreter','latex'); xlabel(xLabel,'interpreter','latex'); set(gca,'fontsize',12); grid on
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')


    fig = figure(5); % Magnetic Field Components and Magnitude in ECI
    plot(plotTimeVec, magFieldMatrix_ECI(:,1), 'Color',[0 .5 .5],'linewidth',3);
    hold on
    plot(plotTimeVec, magFieldMatrix_ECI(:,2), 'Color',[.5 0 .5],'linewidth',3);
    plot(plotTimeVec, magFieldMatrix_ECI(:,3), 'Color',[0.7 0.6 0.3],'linewidth',3);
    plot(plotTimeVec, vecnorm(magFieldMatrix_ECI,2,2), 'Color',[0.1 0.1 0.1],'linewidth',3);
    legend('$ \stackrel{\rightharpoonup}{b}|_{I,x} $', '$ \stackrel{\rightharpoonup}{b}|_{I,y} $', '$ \stackrel{\rightharpoonup}{b}|_{I,z} $', '$ |\stackrel{\rightharpoonup}{b}| $','interpreter','latex','Location','best','FontSize',16)
    titleName ='Magnetic Field on Orbit Resolved in ECI';
    %     title(titleName);
    ylabel('Magnetic Flux Density (nT)','interpreter','latex'); xlabel(xLabel,'interpreter','latex'); set(gca,'fontsize',12); grid on
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')

    fig = figure(6); % Magnetic Field Components and Magnitude in Perifocal Frame
    plot(plotTimeVec, magFieldMatrix_P(:,1), 'Color',[0 .5 .5],'linewidth',3);
    hold on
    plot(plotTimeVec, magFieldMatrix_P(:,2), 'Color',[.5 0 .5],'linewidth',3);
    plot(plotTimeVec, magFieldMatrix_P(:,3), 'Color',[0.7 0.6 0.3],'linewidth',3);
    plot(plotTimeVec, vecnorm(magFieldMatrix_P,2,2), 'Color',[0.1 0.1 0.1],'linewidth',3);
    legend('$ \stackrel{\rightharpoonup}{b}|_{P,x} $', '$ \stackrel{\rightharpoonup}{b}|_{P,y} $', '$ \stackrel{\rightharpoonup}{b}|_{P,z} $', '$ |\stackrel{\rightharpoonup}{b}| $','interpreter','latex','Location','best','FontSize',16)
    titleName ='Magnetic Field on Orbit Resolved in P';
    %     title(titleName);
    ylabel('Magnetic Flux Density (nT)','interpreter','latex'); xlabel(xLabel,'interpreter','latex'); set(gca,'fontsize',12); grid on
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')

    % 3D Magnetic Field

    fig = figure(7);
    magFieldMatrixUnit_NED = magFieldMatrix_NED./vecnorm(magFieldMatrix_NED,2,2);
    quiver3(Xorb(1:100:end,1), Xorb(1:100:end,2), Xorb(1:100:end,3),...
        magFieldMatrixUnit_NED(1:100:end,1),magFieldMatrixUnit_NED(1:100:end,2),magFieldMatrixUnit_NED(1:100:end,3),...
        'Color',[0.1 0.1 0.1],'linewidth',2)

    [XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
    hold on;
    surf(XS*Re, YS*Re, ZS*Re);
    axis equal
    xlabel('x (km)','interpreter','latex');ylabel('y (km)','interpreter','latex');zlabel('z (km)','interpreter','latex'); set(gca,'fontsize',12)
    titleName = 'Magnetic Field on Orbit in 3D Resolved in NED';
    %     title(titleName)
    legend('$ \stackrel{\rightharpoonup}{b}$','interpreter','latex','Location','best','FontSize',16)
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    hold off
    saveas(fig,titleName,'pdf')

    fig = figure(8);
    magFieldMatrixUnit_ECEF = magFieldMatrix_ECEF./vecnorm(magFieldMatrix_ECEF,2,2);
    quiver3(Xorb(1:100:end,1), Xorb(1:100:end,2), Xorb(1:100:end,3),...
        magFieldMatrixUnit_ECEF(1:100:end,1),magFieldMatrixUnit_ECEF(1:100:end,2),magFieldMatrixUnit_ECEF(1:100:end,3),...
        'Color',[0.1 0.1 0.1],'linewidth',2)
    [XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
    hold on;
    surf(XS*Re, YS*Re, ZS*Re);
    axis equal
    xlabel('x (km)','interpreter','latex');ylabel('y (km)','interpreter','latex');zlabel('z (km)','interpreter','latex'); set(gca,'fontsize',12)
    titleName = 'Magnetic Field on Orbit in 3D Resolved in ECEF';
    %     title(titleName)
    legend('$ \stackrel{\rightharpoonup}{b}$','interpreter','latex','Location','best','FontSize',16)
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    hold off
    saveas(fig,titleName,'pdf')

    fig = figure(9);
    magFieldMatrixUnit_ECI = magFieldMatrix_ECI./vecnorm(magFieldMatrix_ECI,2,2);
    quiver3(Xorb(1:100:end,1), Xorb(1:100:end,2), Xorb(1:100:end,3),...
        magFieldMatrixUnit_ECI(1:100:end,1),magFieldMatrixUnit_ECI(1:100:end,2),magFieldMatrixUnit_ECI(1:100:end,3),...
        'Color',[0.1 0.1 0.1],'linewidth',2)
    [XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
    hold on;
    surf(XS*Re, YS*Re, ZS*Re);
    axis equal
    xlabel('x (km)','interpreter','latex');ylabel('y (km)','interpreter','latex');zlabel('z (km)','interpreter','latex'); set(gca,'fontsize',12)
    titleName = 'Magnetic Field on Orbit in 3D Resolved in ECI';
    %     title(titleName)
    legend('$ \stackrel{\rightharpoonup}{b}$','interpreter','latex','Location','best','FontSize',16)
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    hold off
    saveas(fig,titleName,'pdf')

end
if flagPlotDipMoment == 1

    % plot magnetic dipole moments
    fig = figure(10);
    subplot(3,1,1)
    plot(plotTimeVec,dipMomMatrix(:,1),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ m_{c,x} (Am^2)$','interpreter','latex'); grid on;

    subplot(3,1,2)
    plot(plotTimeVec,dipMomMatrix(:,2),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ m_{c,y} (Am^2)$','interpreter','latex'); grid on;

    subplot(3,1,3)
    plot(plotTimeVec,dipMomMatrix(:,3),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ m_{c,z} (Am^2)$','interpreter','latex'); grid on;

    hold off
    if FlagControl == 0
        titleName = 'Magnetic Dipole Moment No Control';
    elseif FlagControl == 1
        titleName = 'Magnetic Dipole Moment Saturated Control';
    elseif FlagControl == 2
        titleName = 'Magnetic Dipole Moment Unsaturated Control';
    end
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')

end

AllEuler321 = [];

if FlagAxisPlot ==1
    AllDCM = [];
end


for k=1:1:length(timeVec)
    DCM = qToDCM(Xorb(k,10:13));
    if FlagAxisPlot ==1
        AllDCM = cat(3,AllDCM,DCM);
    end
    Euler321 = DCM2Angle_321(DCM);
    AllEuler321(:,k) = Euler321;
end

AllEuler321 = AllEuler321';

if FlagEulerPlot == 1
    % plot euler angles
    fig = figure(11);
    subplot(3,1,1)
    plot(plotTimeVec,AllEuler321(:,1).*180/pi,'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \phi (deg)$','interpreter','latex'); grid on;

    subplot(3,1,2)
    plot(plotTimeVec,AllEuler321(:,2).*180/pi,'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \theta (deg)$','interpreter','latex'); grid on;

    subplot(3,1,3)
    plot(plotTimeVec,AllEuler321(:,3).*180/pi,'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ \psi (deg)$','interpreter','latex'); grid on;
    hold off

    if FlagControl == 0
        titleName = 'SC Body Euler Angles No Control';
    elseif FlagControl == 1
        titleName = 'SC Body Euler Angles Saturated Control';
    elseif FlagControl == 2
        titleName = 'SC Body Euler Angles Unsaturated Control';
    end
    %     sgtitle(titleName)
    %         set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')
end

if FlagAxisPlot == 1
    fig = figure(12);
    for ii = floor(linspace(1,length(timeVec)-length(timeVec)/4,4))
        DCMplot = squeeze(AllDCM(:,:,ii));
        i = DCMplot*[1 0 0]'; % pembe
        j = DCMplot*[0 1 0]'; % turkuaz
        k = DCMplot*[0 0 1]'; % yesil
        quiver3(Xorb(ii,1),Xorb(ii,2),Xorb(ii,3),i(1)*6000,i(2)*6000,i(3)*6000,'Color',[0.6 0.4 0.6],'linewidth',2)
        hold on
        quiver3(Xorb(ii,1),Xorb(ii,2),Xorb(ii,3),j(1)*6000,j(2)*6000,j(3)*6000,'Color',[0.03 0.5 0.5],'linewidth',2)
        quiver3(Xorb(ii,1),Xorb(ii,2),Xorb(ii,3),k(1)*6000,k(2)*6000,k(3)*6000,'Color',[0.8    0.7    0.09],'linewidth',2)
    end
    grid on
    plot3(Xorb(:,1), Xorb(:,2), Xorb(:,3),'k-','linewidth',3);
    [XS, YS, ZS] = sphere(30); % plot the Earth
    surf(XS*Re, YS*Re, ZS*Re,'FaceAlpha',0.2,'edgecolor','none');
    xlabel('x (km)','interpreter','latex');ylabel('y (km)','interpreter','latex');zlabel('z (km)','interpreter','latex'); set(gca,'fontsize',12)
    titleName = 'Random sc body frames';
    axis equal
    legend('$ i_B $','$ j_B $','$ k_B $','interpreter','latex','Location','best','FontSize',16)
    %     set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    hold off
    saveas(fig,titleName,'pdf')

end

if flagQuaternionPlot == 1
    q = Xorb(:,10:13);         % qaternion:[q_v q_4]

    % plot quaternions
    fig = figure(13);
    subplot(4,1,1)
    plot(plotTimeVec,q(:,1),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ q_1 $','interpreter','latex'); grid on;

    subplot(4,1,2)
    plot(plotTimeVec,q(:,2),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ q_2 $','interpreter','latex'); grid on;

    subplot(4,1,3)
    plot(plotTimeVec,q(:,3),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ q_3 $','interpreter','latex'); grid on;

    subplot(4,1,4)
    plot(plotTimeVec,q(:,4),'LineWidth',2,'Color',[0 .5 .5])
    xlabel(xLabel,'interpreter','latex'); ylabel('$ q_4 $','interpreter','latex'); grid on;
    hold off
    if FlagControl == 0
        titleName = 'SC Quaternion No Control';
    elseif FlagControl == 1
        titleName = 'SC Quaternion Saturated Control';
    elseif FlagControl == 2
        titleName = 'SC Quaternion Unsaturated Control';
    end
    %     set(gcf,'PaperOrientation','landscape')
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    saveas(fig,titleName,'pdf')

end


% % torques
% for ii=1:length(timeVec)
%         [~,~,T] = scDyn13(timeVec(ii),Xorb(ii,:)',TimeStamp0,FlagMagField,FlagControl,constRot);
%         TT(:,ii) = T;
% end