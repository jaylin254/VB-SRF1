% J.M.C.Clark. Shifted Rayleigh filter: A new algorithm for bearing-only tracking
% 仿真scenario 1

function [X_k,P_k,lossflag]=SRF1_clutter(X0,P0,z_k1,R_k,utm)
global h B pc lossnum
% shifted rayleigh filter
% model : X_k=Fai*X_k_1+V_k
% X=(x,vx)
% z_k1 is the conditional angle ,like bearing
% R_k is the variance of z_k1

Vol=B;
Fai=[1 h;
    0 1];
G=[h*h/2 h]';
var_xt=0.01;

Qs=G*G'*var_xt;
Qtr=eye(2);
H=[1 0;
    0 0];
uts=[0 0]';

% projected measurement 
bk0=[cos(z_k1),sin(z_k1)]';
mk=size(bk0,2); 
% pc=1-1/mk;
% pc=0.5;

X_k_1=X0;
P_k_1 = P0; % initial estimation covariance

%filter

%prediction
X_k_k_1=Fai*X_k_1+uts;
P_k_k_1=Fai*P_k_1*Fai'+Qs;

%covariance of projected measurement bk
Qtm=Qtr+(norm(H*X_k_k_1+utm,2).^2+trace(H*P_k_k_1*H'))*R_k*eye(2);
S_k=H*P_k_k_1*H'+Qtm;

 % %%%%%%%%%%%%%%%%%%%%%%% compute weights %%%%%%%%%%%%%%%%%%%%%%
     
y_k_hat=H*X_k_k_1+utm;    

n=mk;
bk=bk0;
    
for j=1:n
    if isempty(find(eig(S_k)<=0))==0  % 非正定
        lossnum=lossnum+1;
        X_k=X_k_k_1;
        P_k=P_k_k_1;
        lossflag=1;
        break;
    else
        
        S_k_sqrt=chol(S_k,'lower'); % S_k_sqrt'*S_k_sqrt=S_k
        S_k_sqrt_inv=inv(S_k_sqrt);
    %     S_k_sqrt_inv=chol(inv(S_k),'lower');
        s11= S_k_sqrt_inv(1,1);
        s12= S_k_sqrt_inv(1,2);
        s21= S_k_sqrt_inv(2,1);
        s22= S_k_sqrt_inv(2,2);


        g_theta=1/sqrt(bk(:,j)'/S_k*bk(:,j));
        h_theta=bearing_generate(s21*bk(1,j)+s22*bk(2,j),s11*bk(1,j)+s12*bk(2,j),0);
        z=y_k_hat'*S_k_sqrt_inv*[cos(h_theta) sin(h_theta)]';
        f_theta=g_theta*g_theta/sqrt(det(S_k))*1/2/pi*exp(-0.5*(y_k_hat'/S_k*y_k_hat-z*z))*(exp(-0.5*z*z)+sqrt(2*pi)*z*0.5*erfc(-z/sqrt(2)));

        c=1/((1-pc)*f_theta+pc/Vol); %normalizing factor
        q0=c*(1-pc)*f_theta;
        q1=c*pc/Vol;

        % %%%%%%%%%%%%%%%%%%%%  compute X_k_0 %%%%%%%%%%%%%%%%%

        %correction steps
        K_k=P_k_k_1*H'/S_k;
        epsilon_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*bk(:,j)'/S_k*(H*X_k_k_1+utm);

        %compute rou2(epsilon_k)
        Fnormal=0.5*erfc(-1*epsilon_k/sqrt(2));
        rou2=(epsilon_k*exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k*epsilon_k+1)*Fnormal)/(exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k)*Fnormal);

        gamma_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*rou2;
        delta_k=1/(bk(:,j)'/S_k*bk(:,j))*(2+epsilon_k*rou2-rou2*rou2);

        %estimation
        X_k_0=(eye(2)-K_k*H)*X_k_k_1-K_k*utm+gamma_k*K_k*bk(:,j);
        P_k_0=(eye(2)-K_k*H)*P_k_k_1+delta_k*K_k*bk(:,j)*bk(:,j)'*K_k';
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % fusion

        X_k=q0*X_k_0+q1*X_k_k_1;
        P_k=q0*(P_k_0+(X_k-X_k_0)*(X_k-X_k_0)')+q1*(P_k_k_1+(X_k-X_k_k_1)*(X_k-X_k_k_1)');

        X_k_k_1=X_k;
        P_k_k_1=P_k;
        lossflag=0;
    end 
end