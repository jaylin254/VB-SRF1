% J.M.C.Clark. Fusion 2005
% 仿真scenario 
% target state with 12 dimention
% using PDA-SRCKF

function [X_k,P_k,lossflag]=PDA_EKF12_clutter(X0,P0,Z_k,utm)
global h q1 q2 q3 sigma sigma_prime gamma lossnum
global B
% shifted rayleigh filter
% model : X_k=Fai*X_k_1+V_k
% X=12 dimentional
% z_k1 is the conditional angle ,like bearing
% R_k is the variance of z_k1
Pd=1;                                                   %检测概率，当不取1时，后面的a计算出来都是0
Pg=0.9997;  
Vol=B;
lamda=4/Vol;

D=[zeros(2,4) eye(2) eye(2) eye(2) zeros(2)]';
A=[1 h;
    0 1];
Fai=blkdiag(A,A,eye(2),eye(2),eye(2),eye(2))+[zeros(12,10) D]*h;

R=[h*h*h/3 h*h/2;
   h*h/2   h];
Qs=blkdiag(q1*R,q1*R,q2*eye(2),q2*eye(2),q2*eye(2),q3*eye(2));

X_k_1=X0;
P_k_1 = P0; % initial estimation covariance

%filter

%prediction
X_k_k_1=Fai*X_k_1;
P_k_k_1=Fai*P_k_1*Fai'+Qs;


% update
for j=1:3
    
    [H,hh]=Jacobi(X_k_k_1,utm([2*j-1;2*j]),j);
    S_k=H*P_k_k_1*H'+sigma*sigma;
    K_k=P_k_k_1*H'/S_k;
    X_k=X_k_k_1+K_k*(Z_k{j}-hh); 
    P_k=P_k_k_1-K_k*S_k*K_k';
    
    X_k_k_1=X_k;
    P_k_k_1=P_k;   
      
end

for j=4:6
    
    [H,hh]=Jacobi(X_k_k_1,0,j);
    S_k=H*P_k_k_1*H'+sigma_prime*sigma_prime;
    K_k=P_k_k_1*H'/S_k;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    nc=size(Z_k{j},1)-1;
    m=1;
    for i=1:nc+1
        v(i)= Z_k{j}(i)-hh;
        d_squa(i)=v(i)'/S_k*v(i);
        if d_squa(i)<=gamma
           gate_meas(m)=Z_k{j}(i);
           m=m+1;
        end
    end
    if m==1  % no measurement falls in the gate
       X_k=X_k_k_1;% using the predicted state as the renewed state
       P_k=P_k_k_1;
    else
        nc=size(gate_meas,2)-1; % the number of measurements

        bb=lamda*sqrt(2*pi*det(S_k))*(1-Pd*Pg)/Pd;
        for m=1:1:nc+1                                         %关联概率的计算
             vgate(m)= gate_meas(m)-hh;
             e(m)=exp(-0.5*vgate(m)'/(S_k)*vgate(m));     
        end

        beta=e./(bb+sum(e));
        beta0=bb/(bb+sum(e));
        
          % 更新值
          PP=0;vv=0;
          for m=1:nc+1
              vv=vv+beta(m).*vgate(m);
              PP=PP+beta(m).*vgate(m)*vgate(m)';
          end
          X_k=X_k_k_1+K_k*vv;

          P_tilt=K_k*(PP-vv*vv')*K_k';
          P_k=beta0*P_k_k_1+(1-beta0)*(eye(12)-K_k*H)*P_k_k_1+P_tilt;
          
          lossflag=0;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_k_k_1=X_k;
    P_k_k_1=P_k;       
end
   

    