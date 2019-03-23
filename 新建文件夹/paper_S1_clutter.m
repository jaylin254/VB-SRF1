% J.M.C.Clark. Shifted Rayleigh filter: A new algorithm for bearing-only tracking
% 仿真scenario 1

clc;clear; close all;
global h  lossnum lossnum_VB
h=1; % samPDASRF period
F=[1 h;
    0 1];
G=[h*h/2 h]';
var_xt=0.01;
var_xp=1; % variance of the platform position noise
var_eta=0.05*0.05; % variance of sensor angle noise
lossnum=0;lossnum_VB=0;
%%%%%%%%%%%%%%%% parameters for PDASRF
global Pd Pg lamda B gamma pc
Pd=1;
Pg=0.997;
mc=4;%number of clutters
B=360*pi/180; % filed of sight

gamma=9;
Vol=B; % volum of the gate
pc=0.7;
lamda=-log(1-pc)/Vol;
% lamda=mc/B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs=G*G'*var_xt;
Qtr=eye(2);
H=[1 0;
    0 0];
uts=[0 0]';


Xt=[80 1]'; % initial state of the target
%  target track
for k=1:20
    Xt(:,k+1)=F*Xt(:,k);%+sqrt(var_xt)*randn.*G;
end

% platform/sensor track
for k=1:21
    xtp1(k)=4*(k-1)+sqrt(var_xp)*randn;
    xtp2(k)=20+sqrt(var_xp)*randn;
end

figure(1)
yt=zeros(1,21);
plot(Xt(1,:),yt,'r--','LineWidth',2);
hold on
plot(xtp1,xtp2,'b-','LineWidth',2);
hold on
plot(20*ones(1,105),'k:','LineWidth',0.5)
xlabel('Horizontal Displacement','FontSize',12);
ylabel('Vertical Displacement','FontSize',12);
legend('Target Motion','Platform Motion','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小



% run
runs=1000;
nolossnum=1;
for i=1:runs  % monte carlo runs
    
    etac=0.5;
    alpha1=2;
    alpha2=10;
    
    % measurements in angle
    for k=1:21
        
      % 量测产生方法1：第一个量测为真实量测，后边mc个杂波
        eta_m(1,k)=bearing_generate(xtp2(k),Xt(1,k)-xtp1(k),sqrt(var_eta));
        eta_m(2:mc+1,k)=unifrnd(-0.5*B,0.5*B,mc,1);
        
    end
    
    % initialization
    X0_0=[20/tan(eta_m(1)) 0]';
    P0_0=[var_xp*(1+1/tan(eta_m(1)))+400*var_eta/(sin(eta_m(1))^4) 0;
      0                                                        1];
    X_SRF=X0_0;
    P_SRF=P0_0;  
    X_PDASRF=X0_0;
    P_PDASRF=P0_0; 
    X_VBSRF=X0_0;
    P_VBSRF=P0_0; 
    X_PDA=X0_0;
    P_PDA=P0_0; 
    S_PDA=chol(P_PDA,'lower');
    X_MEFPDA=X0_0;
    S_MEFPDA= S_PDA; 
     
    lossflag=0;
    
    for k=1:20
        utm=[-4*k 20]';
%         tic
        [X_VBSRF(:,k+1),P_VBSRF(:,:,k+1),etac(k+1),alpha1(k+1),alpha2(k+1),lossflag]=VB_SRF1_clutter(X_VBSRF(:,k),P_VBSRF(:,:,k),eta_m(:,k+1),var_eta,utm,etac(k),alpha1(k),alpha2(k)); % treatment in paper Fusion 2005
%         toc
        if lossflag==1
           break;
        end 
%         tic
       [X_SRF(:,k+1),P_SRF(:,:,k+1),lossflag]=SRF1_clutter(X_SRF(:,k),P_SRF(:,:,k),eta_m(:,k+1),var_eta,utm);
%        toc
       if lossflag==1
           break;
        end    
%        
%         tic
        [X_PDA(:,k+1),S_PDA(:,:,k+1)]=PDA_SRCKF1_clutter(X_PDA(:,k),S_PDA(:,:,k),eta_m(:,k+1),var_eta,utm);
%         toc
        
%         tic
        [X_MEFPDA(:,k+1),S_MEFPDA(:,:,k+1)]=MEFPDA_SRCKF1_clutter(X_MEFPDA(:,k),S_MEFPDA(:,:,k),eta_m(:,k+1),var_eta,utm);
%         toc
        
       if lossflag==1
           break;
       end
%        
        
    end
    
    if lossflag~=1
        
        r_SRF(nolossnum,:)=(X_SRF(1,:)-Xt(1,:)).^2;
        v_SRF(nolossnum,:)=(X_SRF(2,:)-Xt(2,:)).^2;
        r_VBSRF(nolossnum,:)=(X_VBSRF(1,:)-Xt(1,:)).^2;
        v_VBSRF(nolossnum,:)=(X_VBSRF(2,:)-Xt(2,:)).^2;
        r_PDA(nolossnum,:)=(X_PDA(1,:)-Xt(1,:)).^2;
        v_PDA(nolossnum,:)=(X_PDA(2,:)-Xt(2,:)).^2;
        r_MEFPDA(nolossnum,:)=(X_MEFPDA(1,:)-Xt(1,:)).^2;
        v_MEFPDA(nolossnum,:)=(X_MEFPDA(2,:)-Xt(2,:)).^2;
        
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
plot(eta_m(1,:)/pi*180,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12);
ylabel('Bearing measurement （ ^o ）','FontSize',12);
set(gca,'FontSize',12); % 调整坐标的字体大小
