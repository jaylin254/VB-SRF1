% J.M.C.Clark. Fusion 2005
% 仿真scenario 
% target state with 12 dimention
% using PDA-SRCKF

function [X_k,S_k,lossflag]=MEFPDA_SRCKF12_clutter(X0,S0,Z_k,utm)
global h q1 q2 q3 sigma sigma_prime gamma lamda lossnum lossnum_MEFPDA 
global B 
% shifted rayleigh filter
% model : X_k=Fai*X_k_1+V_k
% X=12 dimentional
% z_k1 is the conditional angle ,like bearing
% R_k is the variance of z_k1
Pd=1;                                                   %检测概率，当不取1时，后面的a计算出来都是0
Pg=0.9997;  
Vol=B;
% lamda=3/Vol;
lossflag=0;

D=[zeros(2,4) eye(2) eye(2) eye(2) zeros(2)]';
A=[1 h;
    0 1];
Fai=blkdiag(A,A,eye(2),eye(2),eye(2),eye(2))+[zeros(12,10) D]*h;

R=[h*h*h/3 h*h/2;
   h*h/2   h];
Qs=blkdiag(q1*R,q1*R,q2*eye(2),q2*eye(2),q2*eye(2),q3*eye(2));

DimX=12;
epsilon=sqrt(DimX)*[eye(DimX) -eye(DimX)];

X_k_1=X0;
S_k_1 = S0; % initial estimation covariance

%filter

%prediction
% X_k_k_1=Fai*X_k_1;
% P_k_k_1=Fai*S_k_1*S_k_1'*Fai'+Qs;
% S_k_k_1=chol(P_k_k_1,'lower');

%prediction
for i=1:2*DimX
    X_k_1_point(:,i)=S_k_1*epsilon(:,i)+X_k_1;
end
for i=1:2*DimX
    X_k_k_1_point(:,i)=Fai*X_k_1_point(:,i);
end
X_k_k_1=sum(X_k_k_1_point,2)/(2*DimX);

e_k_k_1=(X_k_k_1_point-X_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
S_Q_k_1=chol(Qs,'lower');
[~,S_k_k_1]=qr([e_k_k_1 S_Q_k_1]',0);
S_k_k_1=S_k_k_1';
P_k_k_1=S_k_k_1*S_k_k_1';

% update
for j=1:3
    
    for i=1:2*DimX
        X_k_point(:,i)=S_k_k_1*epsilon(:,i)+X_k_k_1;
    end
    
    for i=1:2*DimX
         Z_k_k_1_point(:,i)=bearing_generate(utm(2*j-1)-X_k_point(4+2*j-1,i),utm(2*j)-X_k_point(4+2*j,i),0);
    end
    Z_k_k_1=sum(Z_k_k_1_point,2)/(2*DimX);

    e_zz_k_k_1=(Z_k_k_1_point-Z_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
    [~,S_zz_k_k_1]=qr([e_zz_k_k_1 sigma]',0);
    S_zz_k_k_1=S_zz_k_k_1';

    e_x_k_k_1=( X_k_point-X_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
   
    P_xz=e_x_k_k_1*e_zz_k_k_1';
    K_k=P_xz/S_zz_k_k_1'/S_zz_k_k_1;
    X_k=X_k_k_1+K_k*(Z_k{j}-Z_k_k_1); 

    [~, S_k]=qr([e_x_k_k_1-K_k*e_zz_k_k_1 K_k*sigma]',0);
    S_k=S_k';
    
    X_k_k_1=X_k;
    S_k_k_1=S_k;
    
end

for j=4:6
    
    for i=1:2*DimX
        X_k_point(:,i)=S_k_k_1*epsilon(:,i)+X_k_k_1;
    end
    
    for i=1:2*DimX
         Z_k_k_1_point(:,i)=bearing_generate((X_k_point(1,i)-X_k_point(2*j-3,i)),(X_k_point(3,i)-X_k_point(2*j-2,i)),0);
    end
    Z_k_k_1=sum(Z_k_k_1_point,2)/(2*DimX);

    e_zz_k_k_1=(Z_k_k_1_point-Z_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
    [~,S_zz_k_k_1]=qr([e_zz_k_k_1 sigma_prime]',0);
    S_zz_k_k_1=S_zz_k_k_1';

    e_x_k_k_1=( X_k_point-X_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
   
    P_xz=e_x_k_k_1*e_zz_k_k_1';
    K_k=P_xz/S_zz_k_k_1'/S_zz_k_k_1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    nc=size(Z_k{j},1)-1;
    m=1;
    for i=1:nc+1
        v(i)= Z_k{j}(i)-Z_k_k_1;
        d_squa(i)=v(i)'/(S_zz_k_k_1*S_zz_k_k_1')*v(i);
        if d_squa(i)<=100000
           gate_meas(m)=Z_k{j}(i);
           m=m+1;
        end
    end
    if m==1  % no measurement falls in the gate
       X_k=X_k_k_1;% using the predicted state as the renewed state
       S_k=S_k_k_1;
    else
        nc=size(gate_meas,2)-1; % the number of measurements

        bb=lamda*sqrt(2*pi*det(S_zz_k_k_1*S_zz_k_k_1'))*(1-Pd*Pg)/Pd;
         %%%%% %关联概率的计算 using maximum entropy fuzzy clustering principle
        for m=1:1:nc+1                                         %关联概率的计算
             vgate(m)= gate_meas(m)-Z_k_k_1;
             d(m)=vgate(m)'/(S_zz_k_k_1*S_zz_k_k_1')*vgate(m);     
        end
        
        eta=0.8;
        alpha=eta/lamda/min(d);

         for m=1:1:nc+1                                        
             ed(m)=exp(-alpha*d(m));  
         end

        beta=ed./(sum(ed));
        beta0=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % 更新值
          PP=0;vv=0;
          for m=1:nc+1
              vv=vv+beta(m).*vgate(m);
              PP=PP+beta(m).*vgate(m)*vgate(m)';
          end
          X_k=X_k_k_1+K_k*vv;

          P_tilt=K_k*(PP-vv*vv')*K_k';
          P_k=P_k_k_1-(1-beta0)*K_k*(S_zz_k_k_1*S_zz_k_k_1')*K_k'+P_tilt;
          
          if isempty(find(eig(P_k)<=0))==0  % 非正定
                lossnum=lossnum+1;
                lossnum_MEFPDA=lossnum_MEFPDA+1;
                X_k=X_k_k_1;
                S_k=S_k_k_1;
                lossflag=1;
                break;
          end
          S_k=chol(P_k,'lower');
          lossflag=0;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_k_k_1=X_k;
    S_k_k_1=S_k;       
end
   

    