% This script computes the torques required for cycloidal motion of 
% four bar mechanism (Inverse Dynamics) and then feeds that torque 
% back to mechanism (Forward Dynamics) to see whether the expected 
% motion is obtained

% For reference to different symbols used in the code please follow 
% the paper "Lagrangian Dynamic Formulation of a Four-Bar Mechanism 
% with Minimal Coordinates" by Chin Pei Tang and to the wikipedia
% article of four bar mechanism

% Select the ode function appropriately as ode45 can be extremely slow
% or the starting config may be singular

% Update the axis limit in case you change the link lengths

clear variables;
close all;
clc;

%% Global Variables shared between different functions
global th0 thf t0 tf thStart thFinal
global L0 L1 L2 L3 Lc1 Lc2 Lc3
global m1 m2 m3 J1 J2 J3 P1 g
global Config

L0 = 10;  L1 = 9;  L2 = 7;  L3 = 3;
m1 = 0.5; m2 = 0.7; m3 = 0.6;
Lc1 = L1/2; Lc2 = L2/2; Lc3 = L3/2;
I1 = m1*L1^2/12;
I2 = m2*L2^2/12;
I3 = m3*L3^2/12;

g = 9.81;
t0 = 0;     tf = 5;

% Config 1/-1 corresponds to elbow up and elbow up configuration of four bar
Config = 1;

% invSteps are the number of points at whom the inverse kinematics is 
% computed
invSteps = 43;    

%% Deciding the range of motion
mech = mechanismCheck();

if strcmp(mech,'NA')
   return 
elseif (strcmp(mech, 'GCC') || strcmp(mech, 'GCR'))
    th0 = pi/3;     thf = pi/3 + 2*pi;
    thStart = th0;
    thFinal = thf;
elseif strcmp(mech,'GRC') 
    th0 = acos((L0^2 + L1^2 - (L2-L3)^2)/(2*L0*L1));     
    thf = acos((L0^2 + L1^2 - (L2+L3)^2)/(2*L0*L1));
    thStart = th0+degtorad(1);
    thFinal = thf-degtorad(1);
elseif strcmp(mech,'GRR') && Config == 1
    th0 = acos((L0^2 + L1^2 - (L2-L3)^2)/(2*L0*L1));     
    thf = acos((L0^2 + L1^2 - (L2+L3)^2)/(2*L0*L1));
    thStart = th0+degtorad(1);
    thFinal = thf-degtorad(1);
elseif strcmp(mech,'GRR') && Config == -1
    th0 = -acos((L0^2 + L1^2 - (L2+L3)^2)/(2*L0*L1));
    thf = -acos((L0^2 + L1^2 - (L2-L3)^2)/(2*L0*L1));
    thStart = th0+degtorad(1);
    thFinal = thf-degtorad(1);
end

%% Inverse Kinematics and Dynamics
invTimeSeries = linspace(t0,tf,invSteps);

% dataInKin is used for storing the inverse kinematics information
dataInKin = zeros(invSteps,9);

% dataInDynamic is used for storing the inverse Dynamics information
dataInDynamic = zeros(invSteps,1);

J1 = 1/2*(m1*Lc1^2 + I1 + m2*L1^2);
J2 = 1/2*(m2*Lc2^2 + I2);
J3 = 1/2*(m3*Lc3^2 + I3);
P1 = m2*L1*Lc2;

for i = 1:invSteps
    data = cycloidal(invTimeSeries(i));
    th = data(1);   dth = data(2);  ddth = data(3);    
    [alpha, phi] = inKin(th);    
    % Velocities
    S1 = (L1/L2)*(sin(phi-th)/sin(alpha - phi));
    S2 = (L1/L3)*(sin(alpha-th)/sin(alpha - phi));
    dalpha = S1*dth;
    dphi = S2*dth;  
    
    % Acceleration
    dS1_th = (L1/L2) * ( sin(alpha-phi)*cos(phi-th)*(S2 - 1)...
        - sin(phi-th)*cos(alpha-phi)*(S1 - S2))/((sin(alpha-phi))^2);
    dS2_th = (L1/L3) * ( sin(alpha-phi)*cos(alpha-th)*(S1 - 1)...
        - sin(alpha-th)*cos(alpha-phi)*(S1 - S2))/((sin(alpha-phi))^2);
    ddalpha = S1*ddth + dS1_th*dth^2;
    ddphi = S2*ddth + dS2_th*dth^2;     
    
    dataInKin(i,:) = [th, alpha, phi, dth,dalpha,dphi, ddth,ddalpha,ddphi];
    
    % Inverse Dynamics
    C1 = cos(th-alpha);
    dC1_th = -sin(th-alpha);
    dC1_alpha = sin(th-alpha);
    dG_th = (-m1*g*Lc1 - m2*g*L1)*cos(th);
    dG_alpha = -m2*g*Lc2*cos(alpha);
    dG_phi = -m3*g*Lc3*cos(phi);

    term1 = 2*(J1 + J2*S1^2 + J3*S2^2 + P1*C1*S1);
    term2 = 2*J2*dS1_th + 2*J3*S2*dS2_th + P1*(C1*dS1_th+S1*(dC1_th + S1*dC1_alpha));
    dataInDynamic(i) = term1*ddth + term2*dth^2 - dG_th ...
        - S1*dG_alpha-S2*dG_phi;
end

%% Forward Dynamics
y0 = [dataInKin(1,1) dataInKin(1,4)];
refine = 2;

options = odeset();
tstart = t0;
timeForward = tstart;
forwardOutput = y0;

for i = 1:50
    try
        [t,y] = ode23s(@(t,theta) forwardDyanmics(t, theta, invTimeSeries, dataInDynamic),...
                                 [tstart tf], y0, options);
    catch
        mechanismCheck();
        msg = 'Vinay : Constraint limit was voilated. Trying to restart Ode with updated tStart';

        warning(msg)
        tstart = tstart + 0.05;

        if tstart > tf
            break
        end
        
        if exist('t','var')
            nt = length(t);
            y0(1) = y(nt,1);
            y0(2) = -0.2*y(nt,2);
            options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
                'MaxStep',t(nt)-t(1));
        else
            options = odeset(options,'InitialStep',0.001,...
                'MaxStep',0.1/i);
            if i == 50
                clc;
                msg1 = 'Vinay : The required torque was not enough to make the mechanism move';
                msg2 = 'Or stuck in a singular well';
                warning(strcat(msg1,msg2))
            end
        end
        
        continue
    end
    
    nt = length(t);
    timeForward = [timeForward; t(2:nt)];
    forwardOutput = [forwardOutput; y(2:nt,:)];

    y0(1) = y(nt,1);
    y0(2) = -0.2*y(nt,2);

    tstart = t(nt);
    if tstart == tf
        break;
    end
        options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
            'MaxStep',t(nt)-t(1));  
end            
%%            
%Animation of mechanism as expected for cycloidal profile
% strTitle = 'Cycloidal Motion used for computing Inverse Dynamics';
% mechanismAnimation(invTimeSeries, dataInKin(:,1), strTitle, 'CycloidalMotion.avi')

% Animation of mechanism motion computed using forward dynamics
if exist('forwardOutput','var')
    strTitle = 'Motion computed using forward dynamics';
    mechanismAnimation(timeForward, forwardOutput(:,1), strTitle, 'ForwardDynamics.avi')
end
%% Comparing the results of Inverse and Forward
comparisonPlots(invTimeSeries, dataInKin, timeForward, forwardOutput);

%Plot of inverse kinematic and inverse dynamics data for debugging purpose
InKinNDynamicPlotter(invTimeSeries, dataInKin, dataInDynamic);

%% Function used for forward Dynamics
function out = forwardDyanmics(t, theta, timeSeries, dataInDynamic)
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
    dG_th = (-m1*g*Lc1 - m2*g*L1)*cos(th);
    dG_alpha = -m2*g*Lc2*cos(alpha);
    dG_phi = -m3*g*Lc3*cos(phi);
    term1 = 2*(J1 + J2*S1^2 + J3*S2^2 + P1*C1*S1);
    term2 = 2*J2*dS1_th + 2*J3*S2*dS2_th + P1*(C1*dS1_th+S1*(dC1_th + S1*dC1_alpha));
       
	tau = interp1(timeSeries,dataInDynamic,t);
    out = zeros(2,1);   
    out(1) = dth;
    out(2) = 1/term1*(tau + dG_th + S1*dG_alpha + S2*dG_phi - term2*dth^2);
end

function [value, isterminal, direction] = myEvent(~, theta)
    global th0 thf
%     disp([t,theta(1)]);

    val1 = theta(1) - th0;
    val2 = thf - theta(1);
    
    value = [val1;val2];
    isterminal = [1;1];   % Stop the integration
    direction  = [0;0];
end


%% Fucntion for cycloidal profile
function cycloidalProfile = cycloidal(t)
    global thStart thFinal t0 tf
    th = thStart + (thFinal-thStart)/(tf-t0) * ( t-t0 -  (tf-t0)/(2*pi) * sin(2*pi*(t-t0)/(tf-t0)));
    dth = (thFinal-thStart)/(tf-t0)*(1 - cos(2*pi*(t-t0)/(tf-t0)));
    ddth = (thFinal-thStart)/(tf-t0)*(2*pi/(tf-t0)*sin(2*pi*(t-t0)/(tf-t0)));
    cycloidalProfile = [th, dth, ddth];
end

%% Determining the type of mechanism
% Wikipidea four bar
% A crank: can rotate a full 360 degrees
% A rocker: can rotate through a limited range of angles which does not include 0° or 180°
% A 0-rocker: can rotate through a limited range of angles which includes 0° but not 180°
% A ?-rocker: can rotate through a limited range of angles which includes 180° but not 0°
function mech = mechanismCheck()
    global L0 L1 L2 L3
    conditionL0  = lengthCheck(L0, [L1, L2, L3]);
    conditionL1  = lengthCheck(L1, [L0, L2, L3]);
    conditionL2  = lengthCheck(L2, [L0, L1, L3]);
    conditionL3  = lengthCheck(L3, [L0, L1, L2]);
    
    mech = 'NA';
    if conditionL0*conditionL1*conditionL2*conditionL3 == -1
        disp('Four Bar Mechanism is not possible')
    else
        T1 = sign(L0+L2-L1-L3); 
        T2 = sign(L3+L0-L1-L2);
        T3 = sign(L3+L2-L1-L0);

        if T1 == -1 && T2 == -1 && T3 == 1
            disp('Grashof : Crank Crank');
            mech = 'GCC';
        elseif T1 == 1 && T2 == 1 && T3 == 1
            disp('Grashof : Crank Rocker')
            mech = 'GCR';
        elseif T1 == 1 && T2 == -1 && T3 == -1
            disp('Grashof : Rocker Crank')
            disp('Needs some more work')
            mech = 'GRC';
        elseif T1 == -1 && T2 == 1 && T3 == -1
            disp('Grashof : Rocker Rocker')
            disp('Needs some more work')
            mech = 'GRR';
        elseif T1 == -1 && T2 == -1 && T3 == -1
            disp('Non-Grashof : 0-Rocker 0-Rocker')
            disp('Not yet implemented')
        elseif T1 == -1 && T2 == 1 && T3 == 1
            disp('Non-Grashof : pi-Rocker pi-Rocker')
            disp('Not yet implemented')
        elseif T1 == 1 && T2 == -1 && T3 == 1
            disp('Non-Grashof : pi-Rocker 0-Rocker')
            disp('Not yet implemented')
        elseif T1 == 1 && T2 == 1 && T3 == -1
            disp('Non-Grashof : 0-Rocker pi-Rocker')
            disp('Not yet implemented')
        else
            disp('Special case')
            disp('Not yet implemented')
        end
    end
end

%% This function is used to ensure that the mechanism is possible
% for given link lengths
function result = lengthCheck(a,B)
    if a > sum(B)
        result = -1;
    else
        result = 1;
    end
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
        if t > tLast+0.05
            tLast = t;
            th = thSeries(i);
            [alpha, phi] = inKin(th);    

            dataA = [0,0; L1*cos(th), L1*sin(th); L1*cos(th)+ L2*cos(alpha), ...
                L1*sin(th)+L2*sin(alpha); L0+L3*cos(phi), L3*sin(phi); L0,0];
            cla
            h = plot(dataA(:,1), dataA(:,2),'-o','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 .6 .6]);
            text(-13,-13, strcat('Time : ',num2str(timeSeries(i)), ' s'))
            title(strTitle,'Interpreter','latex');
            axis([-5  15 -10 10]);
            axis square;
            movegui(h, 'onscreen');
            drawnow;
            writeVideo(vidObj,getframe(gcf,rect));
        end
    end
    close(vidObj);
end

%% Function for comparing theta and theta_dot plots of Inverse Dynamics and 
% Forward Dynamics
function comparisonPlots(timeSeries, dataInKin, timeForward, forwardOutput)
    figure()
    plot(timeSeries, dataInKin(:,1),'-o','LineWidth',2)
    hold on
    plot(timeForward, forwardOutput(:,1),'-x','LineWidth',2)
	xlabel('Time ($sec$)','Interpreter','latex')
    ylabel('$\theta \: (rad)$','Interpreter','latex')
	title('Comparison $\theta$ values','Interpreter','latex')
	legend('Cycloidal Profile','Forward Dynamiscs Result');
    hold off
    
    figure()
    hold on
    plot(timeSeries, dataInKin(:,4),'-o','LineWidth',2)
    hold on
    plot(timeForward, forwardOutput(:,2),'-x','LineWidth',2)
	xlabel('Time ($s$)','Interpreter','latex')
	ylabel('$\dot{\theta} \: (\frac{rad}{sec})$ ','Interpreter','latex')
	title('Comparison $\dot{\theta}$ values','Interpreter','latex')
	legend('Cycloidal Profile','Forward Dynamiscs Result');
	hold off
end

%% plots of inverse kinematic and inverse dynamic
function InKinNDynamicPlotter(timeSeries, dataInKin, dataInDynamic)
    figure()
    ax1 = subplot(3,1,1);
    plot(timeSeries, dataInKin(:,1), 'LineWidth',2);
	title('$\theta$ and its derivatives for cycloidal profiles','Interpreter','latex')
    ylabel('$\theta \:(rad)$','Interpreter','latex')
    ax2 = subplot(3,1,2);
    plot(timeSeries, dataInKin(:,4), 'LineWidth',2);
    ylabel('$\dot{\theta}\: (\frac{rad}{sec})$','Interpreter','latex')
    ax3 = subplot(3,1,3);
    plot(timeSeries, dataInKin(:,7), 'LineWidth',2);
	xlabel('$sec$', 'Interpreter', 'latex')
    ylabel('$\ddot{\theta}\:(\frac{rad}{sec^2})$','Interpreter','latex')
	linkaxes([ax1,ax2,ax3],'x');

    figure()
    ax1 = subplot(3,1,1);
    plot(timeSeries, dataInKin(:,2), 'LineWidth',2);
	title('$\alpha$ and its derivatives for cycloidal profiles','Interpreter','latex')
    ylabel('$\alpha \: (rad)$','Interpreter','latex')
    ax2 = subplot(3,1,2);
    plot(timeSeries, dataInKin(:,5), 'LineWidth',2);
    ylabel('$\dot{\alpha} \: (\frac{rad}{sec})$','Interpreter','latex')
    ax3 = subplot(3,1,3);
    plot(timeSeries, dataInKin(:,8), 'LineWidth',2);
	xlabel('$sec$', 'Interpreter', 'latex')
    ylabel('$\ddot{\alpha} \: (\frac{rad}{sec^2})$','Interpreter','latex')
	linkaxes([ax1,ax2,ax3],'x');

    figure()
    ax1 = subplot(3,1,1);
    plot(timeSeries, dataInKin(:,3), 'LineWidth',2);
	title('$\phi$ and its derivatives for cycloidal profiles','Interpreter','latex')
    ylabel('$\phi\: (rad)$','Interpreter','latex')
    ax2 = subplot(3,1,2);
	plot(timeSeries, dataInKin(:,6), 'LineWidth',2);
    ylabel('$\dot{\phi}\: (\frac{rad}{sec})$','Interpreter','latex')
    ax3 = subplot(3,1,3);
    plot(timeSeries, dataInKin(:,9), 'LineWidth',2);
	xlabel('$sec$', 'Interpreter', 'latex')
    ylabel('$\ddot{\phi}\: (\frac{rad}{sec^2})$','Interpreter','latex')
	linkaxes([ax1,ax2,ax3],'x');
    
    figure()
    plot(timeSeries, dataInDynamic, 'LineWidth',2);
    title('Torque Computed using Inverse Dynamics','Interpreter','latex')
	xlabel('$sec$', 'Interpreter', 'latex')
    ylabel('$\tau \:(N-m) $','Interpreter','latex')
end
