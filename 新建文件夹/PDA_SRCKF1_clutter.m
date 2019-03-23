function [X_PDA_update,S_PDA_update]=PDA_SRCKF1_clutter(X_k_1,S_k_1,Z_PDA,R_k,utm)
global h B  lossnum lamda
gama=16; 
                                    
Pd=1;                                                   %检测概率，当不取1时，后面的a计算出来都是0
Pg=0.9997;       
Vol=B;
Fai=[1 h;
    0 1];
G=[h*h/2 h]';
var_xt=0.01;

Qs=G*G'*var_xt;
Qtr=eye(2);
uts=[0 0]';

nc=size(Z_PDA,1)-1;

% lamda=3/Vol;  % 杂波密度  

DimX=2;
epsilon=sqrt(DimX)*[eye(DimX) -eye(DimX)];

%prediction
for i=1:2*DimX
    X_k_1_point(:,i)=S_k_1*epsilon(:,i)+X_k_1;
end
for i=1:2*DimX
    X_k_k_1_point(:,i)=Fai*X_k_1_point(:,i)+uts;
end
X_k_k_1=sum(X_k_k_1_point,2)/(2*DimX);

Q_k=0.01*[0.25 0;
          0 1];

e_k_k_1=(X_k_k_1_point-X_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
S_Q_k_1=chol(Q_k,'lower');
[~,S_k_k_1]=qr([e_k_k_1 S_Q_k_1]',0);
S_k_k_1=S_k_k_1';
P_k_k_1=S_k_k_1*S_k_k_1';

% X_k_k_1=Fai*X_k_1+uts;
% P_k_k_1=Fai*P_k_1*Fai'+Qs;
% S_k_k_1=chol(P_k_k_1,'lower');

% measurement update

for i=1:2*DimX
    X_k_point(:,i)=S_k_k_1*epsilon(:,i)+X_k_k_1;
end

for i=1:2*DimX
    [~, Z_k_k_1_point(:,i)]=Jacobi1(X_k_point(:,i),utm);
end
Z_k_k_1=sum(Z_k_k_1_point,2)/(2*DimX);

e_zz_k_k_1=(Z_k_k_1_point-Z_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);
S_R_k=chol(R_k);
[~,S_zz_k_k_1]=qr([e_zz_k_k_1 S_R_k]',0);
S_zz_k_k_1=S_zz_k_k_1';

e_x_k_k_1=( X_k_point-X_k_k_1*ones(1,2*DimX))/sqrt(2*DimX);

% P_z=S_zz_k_k_1*S_zz_k_k_1';
P_xz=e_x_k_k_1*e_zz_k_k_1';
K_k=P_xz/S_zz_k_k_1'/S_zz_k_k_1;

j=1;
for i=1:nc+1
    v(i)= Z_PDA(i)-Z_k_k_1;
    d_squa(i)=v(i)'/(S_zz_k_k_1*S_zz_k_k_1')*v(i);
    if d_squa(i)<=1000000000
       gate_meas(j)=Z_PDA(i);
       j=j+1;
    end
end
if j==1  % no measurement falls in the gate
   X_PDA_update=X_k_k_1;% using the predicted state as the renewed state
   S_PDA_update=S_k_k_1;
else
    nc=size(gate_meas,2)-1; % the number of measurements

    bb=lamda*sqrt(2*pi*det(S_zz_k_k_1*S_zz_k_k_1'))*(1-Pd*Pg)/Pd;
    for j=1:1:nc+1                                         %关联概率的计算
         vgate(j)= gate_meas(j)-Z_k_k_1;
         e(j)=exp(-0.5*vgate(j)'/(S_zz_k_k_1*S_zz_k_k_1')*vgate(j));     
    end

    beta=e./(bb+sum(e));
    beta0=bb/(bb+sum(e));

      % 更新值
      PP=0;vv=0;
      for j=1:nc+1
          vv=vv+beta(j).*vgate(j);
          PP=PP+beta(j).*vgate(j)*vgate(j)';
      end
      X_PDA_update=X_k_k_1+K_k*vv;

      P_tilt=K_k*(PP-vv*vv')*K_k';
      P_PDA_update=P_k_k_1-(1-beta0)*K_k*(S_zz_k_k_1*S_zz_k_k_1')*K_k'+P_tilt;
      S_PDA_update=chol(P_PDA_update,'lower');
%     
end 