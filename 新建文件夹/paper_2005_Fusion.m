% J.M.C.Clark. Fusion 2005
% 仿真scenario 

clc;clear; close all;
global h q1 q2 q3 sigma sigma_prime 
h=1; % sample period
T=100; %tracking period
q1=0.16;
q2=1;
q3=0.02;
sigma=0.014;
sigma_prime=0.283;

%%%%%%%%%%%%%%%% parameters for PDASRF
global Pd Pg pc lamda B gamma lossnum lossnum_SRF lossnum_VB lossnum_PDA lossnum_MEFPDA
Pd=0.9;
Pg=0.997;
mc=2;%number of clutters
B=360*pi/180; % filed of sight
% lamda=mc/B;
% lamda=0.0000005;
% lamda=0.1;
gamma=9;
Vol=B; % volum of the gate
% pc=1-1/(mc+1)*Pd;
% pc=0.667;
% Pd=1-pc^(mc+1);
pc=0.3;% the prob. that the sensor-target bearing is clutter
lamda=-log(1-pc)/Vol;
% lamda=0.8;

lossnum=0;
lossnum_VB=0;
lossnum_SRF=0;
lossnum_PDA=0;
lossnum_MEFPDA=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=[zeros(2,4) eye(2) eye(2) eye(2) zeros(2)]';
A=[1 h;
    0 1];
F=blkdiag(A,A,eye(2),eye(2),eye(2),eye(2))+[zeros(12,10) D]*h;

R=[h*h*h/3 h*h/2;
   h*h/2   h];
Qs=blkdiag(q1*R,q1*R,q2*eye(2),q2*eye(2),q2*eye(2),q3*eye(2));

H(:,:,1)=[0 0 0 0 -1 0 0 0 0 0 0 0;
          0 0 0 0 0 -1 0 0 0 0 0 0];
H(:,:,2)=[ 0 0 0 0 0 0 -1 0 0 0 0 0;
           0 0 0 0 0 0 0 -1 0 0 0 0];
H(:,:,3)=[ 0 0 0 0 0 0 0 0 -1 0 0 0;
           0 0 0 0 0 0 0 0 0 -1 0 0];
H(:,:,4)=[ 1 0 0 0 -1 0 0 0 0 0 0 0;
           0 0 1 0 0 -1 0 0 0 0 0 0];
H(:,:,5)=[ 1 0 0 0 0 0 -1 0 0 0 0 0;
           0 0 1 0 0 0 0 -1 0 0 0 0];
H(:,:,6)=[ 1 0 0 0 0 0 0 0 -1 0 0 0;
           0 0 1 0 0 0 0 0 0 -1 0 0];
       
Xm0=[0 -300]'; % initial state of the monitoring sensor
%  monitoring sensor platform track
for k=1:T
    Xm(1,k)=0;
    Xm(2,k)=Xm0(2)+6*k;
end

% sonobuoys  track

X1t0=[100,-400]';
X2t0=[-200,-400]';
X3t0=[-200,0]';
for k=1:T
     X1t(1,k)=X1t0(1)+2*k+sqrt(q2)*randn;
     X1t(2,k)=X1t0(2)+1*k+sqrt(q2)*randn;
     
     X2t(1,k)=X2t0(1)+2*k+sqrt(q2)*randn;
     X2t(2,k)=X2t0(2)+1*k+sqrt(q2)*randn;
     
     X3t(1,k)=X3t0(1)+2*k+sqrt(q2)*randn;
     X3t(2,k)=X3t0(2)+1*k+sqrt(q2)*randn;
end

% Brownian motion/ Wiener process

dW1 = sqrt(h)*randn(1,T);   % increments
dW2 = sqrt(h)*randn(1,T);   % increments

X1t(1,:) = X1t(1,:)+cumsum(dW1);             % cumulative sum
X1t(2,:)= X1t(2,:)+cumsum(dW2);             % cumulative sum

dW1 = sqrt(h)*randn(1,T);   % increments
dW2 = sqrt(h)*randn(1,T);   % increments

X2t(1,:) = X2t(1,:)+cumsum(dW1);             % cumulative sum
X2t(2,:)= X2t(2,:)+cumsum(dW2);             % cumulative sum

dW1 = sqrt(h)*randn(1,T);   % increments
dW2 = sqrt(h)*randn(1,T);   % increments

X3t(1,:) = X3t(1,:)+cumsum(dW1);             % cumulative sum
X3t(2,:)= X3t(2,:)+cumsum(dW2);             % cumulative sum

%  plot([1:dt:T],[X1t(1,:)],'r-')   % plot W against t


% target track
Ft=[1 h 0 0;
    0 1 0 0;
    0 0 1 h;
    0 0 0 1];
X0t0=[300 -6 -300 4]';
X0t=Ft*X0t0;
for k=2:T
     X0t(:,k)=Ft*X0t(:,k-1);%+sqrt(q1)/h*[[h*h/2;h]*randn;[h*h/2;h]*randn];
     
end



% run
runs=300; 
% initialization
X0_0=[340 -6 -350 6 100 -400 -200 -400 -200 0 2 2]';
% P0_0=diag([10000 1 10000 1 0 0 0 0 0 0 0.02 0.02]);
P0_0=diag([10000 1 10000 1 0.001 0.001 0.001 0.001 0.001 0.001 0.02 0.02]);
X_SRF=X0_0;
P_SRF=P0_0;  
% X_PLE=X0_0;
% P_PLE=P0_0; 
X_VBSRF=X0_0;
P_VBSRF=P0_0; 
X_PDA=X0_0;
P_PDA=diag([10000 1 10000 1 0.001 0.001 0.001 0.001 0.001 0.001 0.02 0.02]);
S_PDA=chol(P_PDA,'lower');
X_MEFPDA=X0_0;
S_MEFPDA=S_PDA;

nolossnum=1;
for i=1:runs  % monte carlo runs
    
    etac=[0.3 0.3 0.3]';
    alpha1=[2 2 2]';
    alpha2=[20 20 20]';
    
    % measurements in angle
    
    for k=1:T
        
        utm=[Xm(:,k);Xm(:,k);Xm(:,k);zeros(6,1)];
        
        theta_1(1,k)=bearing_generate(Xm(1,k)-X1t(1,k),Xm(2,k)-X1t(2,k),sigma);
        theta_2(1,k)=bearing_generate(Xm(1,k)-X2t(1,k),Xm(2,k)-X2t(2,k),sigma);
        theta_3(1,k)=bearing_generate(Xm(1,k)-X3t(1,k),Xm(2,k)-X3t(2,k),sigma);
        eta_1(1,k)=bearing_generate(X0t(1,k)-X1t(1,k),X0t(3,k)-X1t(2,k),0);
        eta_2(1,k)=bearing_generate(X0t(1,k)-X2t(1,k),X0t(3,k)-X2t(2,k),0);
        eta_3(1,k)=bearing_generate(X0t(1,k)-X3t(1,k),X0t(3,k)-X3t(2,k),0);
        
               
         %%%%%%%%%%%%%%%% clutter measurements: more than one measurements are considered
 
         % %%%%第一个量测是真实目标
        eta_1(2:mc+1,k)=unifrnd(-0.5*B,0.5*B,mc,1);
        eta_2(2:mc+1,k)=unifrnd(-0.5*B,0.5*B,mc,1);
        eta_3(2:mc+1,k)=unifrnd(-0.5*B,0.5*B,mc,1);
        p1=binornd(1,Pd);
        eta_1(1,k)=unifrnd(-0.5*B,0.5*B)*(1-p1)+(eta_1(1,k)+randn*sigma_prime)*p1;
        p2=binornd(1,Pd);
        eta_2(1,k)=unifrnd(-0.5*B,0.5*B)*(1-p2)+(eta_2(1,k)+randn*sigma_prime)*p2;
        p3=binornd(1,Pd);
        eta_3(1,k)=unifrnd(-0.5*B,0.5*B)*(1-p3)+(eta_3(1,k)+randn*sigma_prime)*p3;
       
        Z(:,:,k)={theta_1(:,k) theta_2(:,k) theta_3(:,k) eta_1(:,k) eta_2(:,k) eta_3(:,k)};  % cell数组
         
%         tic
       [X_SRF(:,k+1),P_SRF(:,:,k+1),lossflag]=SRF12_clutter(X_SRF(:,k),P_SRF(:,:,k),Z(:,:,k),utm); % treatment in paper Fusion 2005
%         toc
       if lossflag==1
           break;
       end
%         tic
       [X_VBSRF(:,k+1),P_VBSRF(:,:,k+1),etac(:,k+1),alpha1(:,k+1),alpha2(:,k+1),lossflag]=VB_SRF12_clutter(X_VBSRF(:,k),P_VBSRF(:,:,k),Z(:,:,k),utm,etac(:,k),alpha1(:,k),alpha2(:,k)); % treatment in paper Fusion 2005
%        toc
       if lossflag==1
           break;
       end
%        tic
       [X_PDA(:,k+1),S_PDA(:,:,k+1),lossflag]=PDA_SRCKF12_clutter(X_PDA(:,k),S_PDA(:,:,k),Z(:,:,k),utm); % treatemnt as PDA
%        toc
       if lossflag==1
           break;
       end
%        tic
       [X_MEFPDA(:,k+1),S_MEFPDA(:,:,k+1),lossflag]=MEFPDA_SRCKF12_clutter(X_MEFPDA(:,k),S_MEFPDA(:,:,k),Z(:,:,k),utm); % treatemnt as PDA
%        toc
       if lossflag==1
           break;
       end
%         clear eta_1 eta_2 eta_3

%        
    end
    
    if lossflag~=1
        r_SRF(nolossnum,:)=(X_SRF(1,:)-[X0t0(1) X0t(1,:)]).^2+(X_SRF(3,:)-[X0t0(3) X0t(3,:)]).^2;
        v_SRF(nolossnum,:)=(X_SRF(2,:)-[X0t0(2) X0t(2,:)]).^2+(X_SRF(4,:)-[X0t0(4) X0t(4,:)]).^2;
       
        r_VBSRF(nolossnum,:)=(X_VBSRF(1,:)-[X0t0(1) X0t(1,:)]).^2+(X_VBSRF(3,:)-[X0t0(3) X0t(3,:)]).^2;
        v_VBSRF(nolossnum,:)=(X_VBSRF(2,:)-[X0t0(2) X0t(2,:)]).^2+(X_VBSRF(4,:)-[X0t0(4) X0t(4,:)]).^2;
               
        r_PDA(nolossnum,:)=(X_PDA(1,:)-[X0t0(1) X0t(1,:)]).^2+(X_PDA(3,:)-[X0t0(3) X0t(3,:)]).^2;
        v_PDA(nolossnum,:)=(X_PDA(2,:)-[X0t0(2) X0t(2,:)]).^2+(X_PDA(4,:)-[X0t0(4) X0t(4,:)]).^2;
                
        r_MEFPDA(nolossnum,:)=(X_MEFPDA(1,:)-[X0t0(1) X0t(1,:)]).^2+(X_MEFPDA(3,:)-[X0t0(3) X0t(3,:)]).^2;
        v_MEFPDA(nolossnum,:)=(X_MEFPDA(2,:)-[X0t0(2) X0t(2,:)]).^2+(X_MEFPDA(4,:)-[X0t0(4) X0t(4,:)]).^2;

        nolossnum=nolossnum+1;
    end
end
nolossnum=nolossnum-1;

RMSE_r_SRF=sqrt(sum(r_SRF,1)/nolossnum);
RMSE_v_SRF=sqrt(sum(v_SRF,1)/nolossnum);

RMSE_r_VBSRF=sqrt(sum(r_VBSRF,1)/nolossnum);
RMSE_v_VBSRF=sqrt(sum(v_VBSRF,1)/nolossnum);

RMSE_r_PDA=sqrt(sum(r_PDA,1)/nolossnum);
RMSE_v_PDA=sqrt(sum(v_PDA,1)/nolossnum);

RMSE_r_MEFPDA=sqrt(sum(r_MEFPDA,1)/nolossnum);
RMSE_v_MEFPDA=sqrt(sum(v_MEFPDA,1)/nolossnum);


figure(1)
plot(X0t(1,:),X0t(3,:),'k.');
hold on
plot(X_SRF(1,:),X_SRF(3,:),'g--','LineWidth',2);
hold on
plot(X_VBSRF(1,:),X_VBSRF(3,:),'k-','LineWidth',2);
hold on
plot(Xm(1,:),Xm(2,:),'r-.','LineWidth',2);
hold on
plot(X1t(1,:),X1t(2,:),'b-','LineWidth',2);
hold on
plot(X_SRF(5,:),X_SRF(6,:),'g--','LineWidth',2);
hold on
plot(X_VBSRF(5,:),X_VBSRF(6,:),'k-','LineWidth',2);
hold on
legend('True target track','Estimated target track(SRF)','Estimated target track(VBSRF)','Monitoring Sensor','Sensor tracks','Estimated sensor track(SRF)','Estimated sensor track(VBSRF)');
plot(X2t(1,:),X2t(2,:),'b-','LineWidth',2);
hold on
plot(X3t(1,:),X3t(2,:),'b-','LineWidth',2);
hold on
plot(X_SRF(7,:),X_SRF(8,:),'g--','LineWidth',2);
hold on
plot(X_VBSRF(7,:),X_VBSRF(8,:),'k-','LineWidth',2);
hold on
plot(X_SRF(9,:),X_SRF(10,:),'g--','LineWidth',2);
hold on
plot(X_VBSRF(9,:),X_VBSRF(10,:),'k-','LineWidth',2);
xlabel('x position (m)','FontSize',12);
ylabel('y position (m)','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小


figure(2)
plot(RMSE_r_VBSRF,'r','LineWidth',2);
hold on
plot(RMSE_r_SRF,'b-.','LineWidth',2);
hold on
plot(RMSE_r_PDA,'m:','LineWidth',2);
hold on
plot(RMSE_r_MEFPDA,'g--','LineWidth',2);
hold on
xlabel('Time (s)','FontSize',12);
ylabel('RMS position errors (m)','FontSize',12);
legend('VB-SRF','SRF','PDA-SRCKF','MEFPDA-SRCKF','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小


figure(3)
plot(RMSE_v_VBSRF,'r','LineWidth',2);
hold on
plot(RMSE_v_SRF,'b-.','LineWidth',2);
hold on
plot(RMSE_v_PDA,'m:','LineWidth',2);
hold on
plot(RMSE_v_MEFPDA,'g--','LineWidth',2);
hold on
xlabel('Time (s)','FontSize',12);
ylabel('RMS velocity errors (m/s)','FontSize',12);
legend('VB-SRF','SRF','PDA-SRCKF','MEFPDA-SRCKF','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小



figure(4)
plot(eta_1(1,:)*180/pi,'r','LineWidth',2);
hold on
plot(eta_2(1,:)*180/pi,'b--','LineWidth',2);
hold on
plot(eta_3(1,:)*180/pi,'k-.','LineWidth',2);
hold on
xlabel('Time (s)','FontSize',12);
ylabel('Target Measurements ( ^o )','FontSize',12);
legend('eta_1','eta_2','eta_3','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小

