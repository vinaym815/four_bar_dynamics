clear variables;
close all;
clc;

%% Global Variables shared between different functions
global L0 L1 L2 L3 Lc1 Lc2 Lc3
global m1 m2 m3 J1 J2 J3 P1 g
global Config

L0 = 3;  L1 = 9;  L2 = 10;  L3 = 7;
m1 = 1; m2 = 1; m3 = 1;
Lc1 = L1/2; Lc2 = L2/2; Lc3 = L3/2;

t0 = 0;     tf = 15;
Config = -1;
th0 = 1.5708;   dth0 = 0;

g = 9.8;
I1 = m1*L1^2/12;    I2 = m2*L2^2/12;    I3 = m3*L3^2/12;

%% Inverse Kinematics and Dynamics
J1 = 1/2*(m1*Lc1^2 + I1 + m2*L1^2);
J2 = 1/2*(m2*Lc2^2 + I2);
J3 = 1/2*(m3*Lc3^2 + I3);
P1 = m2*L1*Lc2;

%% Forward Dynamics
opt = odeset('RelTol',1e-15,'AbsTol', 1e-15);
[t,y] = ode45(@(t,theta) forwardDyanmics(t, theta),...
                                 [t0 tf], [th0 dth0],opt);

%%   Animation of mechanism motion computed using forward dynamics
strTitle = 'Forward dynamics constant $\tau \; =\; 6\;N-m$';
mechanismAnimation(t, y(:,1), strTitle, 'theValidationVideo.avi')

figure()
plot(t,y(:,1),'LineWidth',3)
title('$\theta \;(rad)$ vs $Time\;(sec)$','Interpreter','latex');
ylabel('$\theta \;(rad)$','Interpreter','latex');
xlabel('$sec$','Interpreter','latex');

hold on
%% Function used for forward Dynamics
function out = forwardDyanmics(~, theta)
    global L1 L2 L3 Lc1 Lc2 Lc3
    global m1 m2 m3 J1 J2 J3 P1 g

    th = theta(1);    dth = theta(2);
    [alpha, phi] = inKin(th);    
        
    S1 = (L1/L2)*(sin(phi-th)/sin(alpha - phi));
    S2 = (L1/L3)*(sin(alpha-th)/sin(alpha - phi));
    
    dS1_th = (L1/L2) * ( sin(alpha-phi)*cos(phi-th)*(S2 - 1)...
        - sin(phi-th)*cos(alpha-phi)*(S1 - S2))/((sin(alpha-phi))^2);
    dS2_th = (L1/L3) * ( sin(alpha-phi)*cos(alpha-th)*(S1 - 1)...
        - sin(alpha-th)*cos(alpha-phi)*(S1 - S2))/((sin(alpha-phi))^2);
    
    C1 = cos(th-alpha);
    dC1_th = -sin(th-alpha);
    dC1_alpha = sin(th-alpha);
    dG_th = (-m1*Lc1 - m2*L1)*g*cos(th);
    dG_alpha = -m2*g*Lc2*cos(alpha);
    dG_phi = -m3*g*Lc3*cos(phi);
    term1 = 2*(J1 + J2*S1^2 + J3*S2^2 + P1*C1*S1);
    term2 = 2*J2*dS1_th + 2*J3*S2*dS2_th + P1*(C1*dS1_th+S1*(dC1_th + S1*dC1_alpha));
    tau = 6;
%     disp([tau,dG_th + S1*dG_alpha + S2*dG_phi - term2*dth^2])
    out = zeros(2,1);   
    out(1) = dth;
    out(2) = (1/term1)*(tau + dG_th + S1*dG_alpha + S2*dG_phi - term2*dth^2);
end


%% InverseKinematics
function [alpha, phi] = inKin(th)
	global L0 L1 L2 L3 Config
    k1 = -2*L1*L3*sin(th);
    k2 = 2*L3*(L0 - L1*cos(th));
    k3 = L0^2 + L1^2 - L2^2 + L3^2 -2*L0*L1*cos(th);
    phi = 2*atan2(-k1 + Config*sqrt(k1^2 + k2^2 - k3^2), k3 - k2);
    alpha = atan2(L3*sin(phi) - L1*sin(th), L0 - L1*cos(th) + L3*cos(phi));
end

%% Animation of mechanism as dictated by forward dynamics
function mechanismAnimation(timeSeries, thSeries, strTitle,fileName)
    global L0 L1 L2 L3 t0
    
    h = figure();
    movegui(h, 'onscreen');
    rect = get(h,'Position'); 
    rect(1:2) = [0 0]; 
    vidObj = VideoWriter(fileName);
    vidObj.Quality = 100;
    open(vidObj);
    tLast = t0-0.12;
    for i = 1:size(thSeries,1)
        t = timeSeries(i);
        if t > tLast + 0.05
            tLast = t;
            th = thSeries(i);
            [alpha, phi] = inKin(th);
%             disp([th,alpha,phi]);

            dataA = [0,0; L1*cos(th), L1*sin(th); L1*cos(th)+ L2*cos(alpha), ...
                L1*sin(th)+L2*sin(alpha); L0+L3*cos(phi), L3*sin(phi); L0,0];
            cla
            h = plot(dataA(:,1), dataA(:,2),'-o','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 .6 .6]);
            text(-13,-13, strcat('Time : ',num2str(timeSeries(i)), ' s'))
            title(strTitle,'Interpreter','latex');
            axis([-15  15 -15 15]);
            axis square;
            movegui(h, 'onscreen');
            drawnow;
            writeVideo(vidObj,getframe(gcf,rect));
        end
    end
    close(vidObj);
end