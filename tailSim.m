function out = tailSim(tail,load,Ts,zeta,ht,ndp,ic,lim)

%-------Inputs-------------%
% tail - 'lcm' for linear counter mass tail
%      - 'arc' for revolute arc tail
% load - numeric value representative of the load being held [kg]
% Ts   - settling time for controller [sec]
% zeta - damping ratio
% ht - height where the tail is mounted on the user's back [m]
% ndp - location of the non-dominant poles
% ic - the initial ankle angle [deg]
% lim  - 0 for no physical limits on tail motion
%      - 1 for limits on tail motion
%
% Example: out = tailSim('lcm',5,4,0.8,0.997,10,8,0) will simulate counter
% mass tail, for a load of 5 kg, with desired controller settling time of 4
% seconds, damping ratio 0.8. Mounted at a height of 0.997. Non dominant
% poles located 10 times away from the dominant poles. At an ic of 8
% degrees. and no limit on the tail motion
%
%-------Outputs------------%
% out  - data structure
%    .lin.staPt - stationary point dynamics linearised about
%    .lin.A     - linearised state matrix
%    .lin.B     - linearised input matrix
%    .CC.A      - control canonical state matrix
%    .CC.B      - control canonical input matrix
%    .gain      - feedback control gains
%    .states    - time solution of state vector
%    .time      - time vector [s]

close all
fig = figure;
ax = axes;

%% Parameters

% Trunk
m0 = 82.2;
h0 = 1.8;
w = 0.15;
mb = m0 - 2*(1.45+0.53);
hb = 0.997; % CoM human standing
I0 = 0*m0*(w^2+h0^2)/12;
%Ib = m0*(w^2+h0^2)/12;
Ib = 12.92;

% Forearm + Hand
ma = 2*(1.45+0.53);
ha = h0-.323-.366; % height - head/neck - shoulder/elbow
la = (0.193+0.257)/2;
wa = .301/pi;
Ia = 0*ma*(wa^2+(2*la)^2)/12;

% Load
ml = load;
dl = .15;
wl = 0.193;
ll = 0.257+.193/2;
hl = ha+wa/2+dl/2;
Il = 0*ml*(wl^2+dl^2)/12;

% Tail
mt = 5;
lt = 0.9;
%ht = hin;
It = mt*lt^2/12;
% Gravity
g = 9.81;

% Ankle Parameters
kp = 10*180/pi; % Stiffness
kd = 4*180/pi;  % Damping

%% Develop State Models
switch tail
    case 'lcm'
        mc = 3;
        Ic = 0;

        % Define stationary point
        th1 = 0;
        dc = (ma*la+ml*ll-mt*lt/2)/mc;
        sp = [th1 dc];

        % Define intermediate variables
        m = mb + ma + ml +mt;
        alpha1 = mb/m;
        alpha2 = ma/m;
        alpha3 = ml/m;
        alpha4 = mt/m;
        
        beta2 = alpha2*la + alpha3*ll - alpha4*lt/2;
        beta1 = alpha1*hb + alpha2*ha + alpha3*hl + alpha4*ht;
        
        m11 = mc*dc^2 + mc*ht^2 + Ib + m*(beta1^2 + beta2^2);
        m12 = mc*ht;
        m22 = mc;
        mu = m11*m22 - m12^2;

        % State Matrices
        A = [-(m22*(kp - g*(mc*(ht*cos(th1) + dc*sin(th1)) + beta1*m*cos(th1) - beta2*m*sin(th1))) - g*m12*mc*cos(th1)), -(g*m22*mc*cos(th1)), -(kd*m22), 0;...
            (m12*(kp - g*(mc*(ht*cos(th1) + dc*sin(th1)) + beta1*m*cos(th1) - beta2*m*sin(th1))) - g*m11*mc*cos(th1)),  (g*m12*mc*cos(th1)),  (kd*m12), 0]/mu;
        Alin = [zeros(2,2) eye(2);A];
        Blin = [0;0;-m12;m11]/mu;
        bound = [0 lt];

    case 'arc'
        % Define stationary point
        th1 = 0;
        th2 = asin((ma*la+ml*ll)/(mt*lt));
        sp = [th1 th2];
        
        % Define intermediate variables
        m = mb + ma + ml;
        alpha1 = mb/m;
        alpha2 = ma/m;
        alpha3 = ml/m;
        
        beta2 = alpha2*la + alpha3*ll;
        beta1 = alpha1*hb + alpha2*ha + alpha3*hl;
        %r = sqrt(beta1^2+beta2^2)
        
        m11 = m*beta1^2 + m*beta2^2 + mt*ht^2 + 2*mt*cos(th2)*ht*lt + mt*lt^2 + Ib + It;
        m12 = (mt*lt^2 + ht*mt*cos(th2)*lt + It);
        m22 = mt*lt^2 + It;
        mu = m11*m22 - m12^2;

        % State Matrices
        A = [(g*lt*m12*mt*cos(th1 - th2)) - (m22*(kp - g*(beta1*m*cos(th1) + ht*mt*cos(th1) - beta2*m*sin(th1) + lt*mt*cos(th1 - th2)))), - (g*lt*m12*mt*cos(th1 - th2)) - (g*lt*m22*mt*cos(th1 - th2)), -(kd*m22), 0;...
             (m12*(kp - g*(beta1*m*cos(th1) + ht*mt*cos(th1) - beta2*m*sin(th1) + lt*mt*cos(th1 - th2)))) - (g*lt*m11*mt*cos(th1 - th2)),   (g*lt*m11*mt*cos(th1 - th2)) + (g*lt*m12*mt*cos(th1 - th2)),  (kd*m12), 0]/mu;
        Alin = [zeros(2,2) eye(2);A];
        Blin = [0;0;-m12;m11]/mu;
        bound = [10 170]*pi/180;
end

% Feedback Control Law
[Acc,Bcc,K] = fullStateFeedback(Alin,Blin,Ts,zeta);

%% Simulation
delT = 0.001;
t = 0:delT:Ts+1;
% initial condition
x0 = [ic*pi/180;0;0;0];
x = x0;
X = x';
for i = 1:length(t)-1
    x_dot = (Alin-Blin*K)*x;
    x = X(i,:)'+x_dot*delT;

    if lim == 1
        if x(2)+sp(2) < bound(1)
            x(2) = bound(1)-sp(2);
        elseif x(2)+sp(2) > bound(2)
            x(2) = bound(2)-sp(2);
        end
    end

    X = [X;x'];
    % Plot
    if mod(i,5)==0
        cla(ax)
        plotPendulum(X(i,1:2)+sp,tail,ax)
    end
end

out.lin.staPt = sp;
out.lin.A = Alin;
out.lin.B = Blin;
out.CC.A = Acc;
out.CC.B = Bcc;
out.gain = K;
out.states = X;
out.time = t;
out.param.h = ht;
out.param.zeta = zeta;
out.param.Ts = Ts;
out.param.ndp = ndp;
out.initial = x0;

%% Nested Functions
    function [Acc,Bcc,K] = fullStateFeedback(A,B,Ts,zeta)

        % Desired Performance
        %ndp = 4;
        wn = 4/(zeta*Ts);
        desCE = conv([1 2*zeta*wn wn^2],[1 2*ndp*zeta*wn (ndp*wn)^2]);
        
        % Convert to Control Canonical Form
        poleOl = eig(A);
        ceOL = 1;
        for ii = 1:length(poleOl)
            ceOL = (conv(ceOL,[1 -poleOl(ii)]));
        end
        ceOL = real(ceOL);
        
        % Control Canonical Form
        Acc = [-ceOL(2:end);...
                1 0 0 0;...
                0 1 0 0;...
                0 0 1 0];
        Bcc = [1;0;0;0];
        
        % Controllability Matrices
        ctrbOrig = [B A*B A^2*B A^3*B];
        ctrbCC = [Bcc Acc*Bcc Acc^2*Bcc Acc^3*Bcc];
        T = ctrbCC/ctrbOrig;

        % Full State Feedback
        Kcc = desCE(2:end)- ceOL(2:end);
        K = Kcc*T;

    end

    function plotPendulum(x,tail,ax)
        %global mb ma ml mt la ll lt h0 ha hb hl ht dl wl
        switch tail
            case 'lcm'
                theta = x(1);
                d = x(2);
                trunkCoM = [hb*sin(theta) hb*cos(theta)];
                trunkCoord = [0 0;h0*sin(theta) h0*cos(theta)];
                plotTrunk = line(ax,trunkCoord(:,1),trunkCoord(:,2));
                plotTrunk.LineWidth = 1.5;
                plotTrunk.Color = [230 81 39]/255;
                
                plotTrunkCoM = line(ax,trunkCoM(:,1),trunkCoM(:,2));
                plotTrunkCoM.LineWidth = 1.5;
                plotTrunkCoM.Marker = 'o';
                plotTrunkCoM.MarkerSize = 10;
                plotTrunkCoM.Color = [230 81 39]/255;
            
                foreArmCoM = [ha*sin(theta)+la*cos(theta) ha*cos(theta)-la*sin(theta)];
                foreArmCoord = [ha*sin(theta) ha*cos(theta);ha*sin(theta)+2*la*cos(theta) ha*cos(theta)-2*la*sin(theta)];
                plotForeArm = line(ax,foreArmCoord(:,1),foreArmCoord(:,2));
                plotForeArm.LineWidth = 1.5;
                plotForeArm.Color = [230 81 39]/255;
                
                plotForeArmCoM = line(ax,foreArmCoM(:,1),foreArmCoM(:,2));
                plotForeArmCoM.LineWidth = 1.5;
                plotForeArmCoM.Marker = 'o';
                plotForeArmCoM.MarkerSize = 10;
                plotForeArmCoM.Color = [230 81 39]/255;
                
                T01 = [cos(pi/2-theta) -sin(pi/2-theta) 0;...
                       sin(pi/2-theta)  cos(pi/2-theta) 0;...
                            0               0       1];
                boxCoord(1,:) = (T01*[ha;-(2*la-wl);1])';
                boxCoord(2,:) = (T01*[ha;-2*la;1])';
                boxCoord(3,:) = (T01*[hl+dl/2;-2*la;1])';
                boxCoord(4,:) = (T01*[hl+dl/2;-(2*la-wl);1])';
                boxCoord(5,:) = boxCoord(1,:);
                plotBox = line(ax,boxCoord(:,1),boxCoord(:,2));
                plotBox.LineWidth = 1.5;
                plotBox.Color = [230 81 39]/255;

                boxCoM   = (T01*[hl;-ll;1])';
                plotBoxCoM = line(ax,boxCoM(1),boxCoM(2));
                plotBoxCoM.LineWidth = 1.5;
                plotBoxCoM.Marker = 'o';
                plotBoxCoM.MarkerSize = 10;
                plotBoxCoM.Color = [230 81 39]/255;
                
                tailCoM = [ht*sin(theta)-lt*cos(theta)/2 ht*cos(theta)+lt*sin(theta)/2];
                tailCoord = [ht*sin(theta) ht*cos(theta);...
                             ht*sin(theta)-lt*cos(theta) ht*cos(theta)+lt*sin(theta)];
                plotTail = line(ax,tailCoord(:,1),tailCoord(:,2));
                plotTail.LineWidth = 1.5;
                plotTail.Color = [230 81 39]/255;
                
                plotTailCoM = line(ax,tailCoM(:,1),tailCoM(:,2));
                plotTailCoM.LineWidth = 1.5;
                plotTailCoM.Marker = 'o';
                plotTailCoM.MarkerSize = 10;
                plotTailCoM.Color = [230 81 39]/255;
                
                counterCoM = [ht*sin(theta)-d*cos(theta) ht*cos(theta)+d*sin(theta)];
                plotCounterCoM = line(ax,counterCoM(:,1),counterCoM(:,2));
                plotCounterCoM.LineWidth = 1.5;
                plotCounterCoM.Marker = 'o';
                plotCounterCoM.MarkerSize = 12;
                plotCounterCoM.Color = [230 81 39]/255;
                
                m = mb+ma+ml+mt;
                alpha_1 = mb/m;
                alpha_2 = ma/m;
                alpha_3 = ml/m;
                alpha_4 = mt/m;
                sysCoM = alpha_1*trunkCoM + alpha_2*foreArmCoM + alpha_3*boxCoM(1:2) + alpha_4*tailCoM;
                
                plotSysCoM = line(ax,sysCoM(:,1),sysCoM(:,2));
                plotSysCoM.LineWidth = 1.5;
                plotSysCoM.Marker = 'o';
                plotSysCoM.MarkerSize = 10;
                plotSysCoM.Color = [100 230 39]/255;
                
                xlim(ax,[-1.01*h0 1.01*h0])
                ylim(ax,[-0.01*h0 1.1*h0])
                
                drawnow
            case 'arc'
                th1p = x(1); th2p = x(2); 

                T01 = [cos(pi/2-th1p) -sin(pi/2-th1p) 0;...
                       sin(pi/2-th1p)  cos(pi/2-th1p) 0;...
                            0               0       1];
                T12 = [cos(th2p)       -sin(th2p)    ht;...
                       sin(th2p)        cos(th2p)    0;...
                            0               0      1];
                T02 = T01*T12;
                
                trunkCoM = T01*[hb;0;1];
                trunkCoord(1,:) = [0 0 1];
                trunkCoord(2,:) = (T01*[h0;0;1])';
                
                plotTrunk = line(ax,trunkCoord(:,1),trunkCoord(:,2));
                plotTrunk.LineWidth = 1.5;
                plotTrunk.Color = [230 81 39]/255;
                
                plotTrunkCoM = line(ax,trunkCoM(1),trunkCoM(2));
                plotTrunkCoM.LineWidth = 1.5;
                plotTrunkCoM.Marker = 'o';
                plotTrunkCoM.MarkerSize = 10;
                plotTrunkCoM.Color = [230 81 39]/255;
                
                foreArmCoord(1,:) = (T01*[ha;0;1])';
                foreArmCoord(2,:) = (T01*[ha;-2*la;1])';
                
                plotForeArm = line(ax,foreArmCoord(:,1),foreArmCoord(:,2));
                plotForeArm.LineWidth = 1.5;
                plotForeArm.Color = [230 81 39]/255;
                
                foreArmCoM = T01*[ha; -la;1];
                plotForeArmCoM = line(ax,foreArmCoM(1),foreArmCoM(2));
                plotForeArmCoM.LineWidth = 1.5;
                plotForeArmCoM.Marker = 'o';
                plotForeArmCoM.MarkerSize = 10;
                plotForeArmCoM.Color = [230 81 39]/255;
                
                boxCoord(1,:) = (T01*[ha;-(2*la-wl);1])';
                boxCoord(2,:) = (T01*[ha;-2*la;1])';
                boxCoord(3,:) = (T01*[hl+dl/2;-2*la;1])';
                boxCoord(4,:) = (T01*[hl+dl/2;-(2*la-wl);1])';
                boxCoord(5,:) = boxCoord(1,:);
                
                plotBox = line(ax,boxCoord(:,1),boxCoord(:,2));
                plotBox.LineWidth = 1.5;
                plotBox.Color = [230 81 39]/255;
                
                boxCoM   = T01*[hl;-ll;1];
                plotBoxCoM = line(ax,boxCoM(1),boxCoM(2));
                plotBoxCoM.LineWidth = 1.5;
                plotBoxCoM.Marker = 'o';
                plotBoxCoM.MarkerSize = 10;
                plotBoxCoM.Color = [230 81 39]/255;
                
                tailCoord(1,:) = T02(:,3)';
                tailCoord(2,:) = (T02*[lt;0;1])';
                plotTail = line(ax,tailCoord(:,1),tailCoord(:,2));
                plotTail.LineWidth = 1.5;
                plotTail.Color = [230 81 39]/255;
                
                tailCoM = T02*[lt;0;1];
                plotTailCoM = line(ax,tailCoM(1),tailCoM(2));
                plotTailCoM.LineWidth = 1.5;
                plotTailCoM.Marker = 'o';
                plotTailCoM.MarkerSize = 10;
                plotTailCoM.Color = [230 81 39]/255;
                
                m = mb+ma+ml;
                alpha_1 = mb/m;
                alpha_2 = ma/m;
                alpha_3 = ml/m;
                
                sysCoM = alpha_1*trunkCoM + alpha_2*foreArmCoM + alpha_3*boxCoM;
                
                plotSysCoM = line(ax,sysCoM(1),sysCoM(2));
                plotSysCoM.LineWidth = 1.5;
                plotSysCoM.Marker = 'o';
                plotSysCoM.MarkerSize = 10;
                plotSysCoM.Color = [100 230 39]/255;
                
                xlim(ax,[-1.01*h0 1.01*h0])
                ylim(ax,[-0.01*h0 1.1*h0])
                
                drawnow
        end
    end
end