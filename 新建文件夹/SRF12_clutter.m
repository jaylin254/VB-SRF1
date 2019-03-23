% J.M.C.Clark. Fusion 2005
% 仿真scenario 
% target state with 12 dimention

function [X_k,P_k,lossflag]=SRF12_clutter(X0,P0,Z_k,utm)
global h q1 q2 q3 sigma sigma_prime pc lossnum lossnum_SRF
global B
% shifted rayleigh filter
% model : X_k=Fai*X_k_1+V_k
% X=12 dimentional
% z_k1 is the conditional angle ,like bearing
% R_k is the variance of z_k1

Vol=B;
D=[zeros(2,4) eye(2) eye(2) eye(2) zeros(2)]';
A=[1 h;
    0 1];
Fai=blkdiag(A,A,eye(2),eye(2),eye(2),eye(2))+[zeros(12,10) D]*h;

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

Qtr=0;
% projected measurement 
% bk=[sin(Z_k(1)),cos(Z_k(1)),sin(Z_k(2)),cos(Z_k(2)),sin(Z_k(3)),cos(Z_k(3)),sin(Z_k(4)),cos(Z_k(4)),sin(Z_k(5)),cos(Z_k(5)),sin(Z_k(6)),cos(Z_k(6))]';


X_k_1=X0;
P_k_1 = P0; % initial estimation covariance

%filter

%prediction
X_k_k_1=Fai*X_k_1;
P_k_k_1=Fai*P_k_1*Fai'+Qs;


%correction step

for i=1:6
        
    bk0=[sin(Z_k{i}),cos(Z_k{i})]';
    mk=size(bk0,2); % number of measurements for one sensor
%     pc=1-1/mk;
   
    if i<=3
        Qtm=Qtr+(norm(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]),2).^2+trace(H(:,:,i)*P_k_k_1*H(:,:,i)'))*sigma*sigma*eye(2);
    else
        Qtm=Qtr+(norm(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]),2).^2+trace(H(:,:,i)*P_k_k_1*H(:,:,i)'))*sigma_prime*sigma_prime*eye(2);
    end
    
       
    % %%%%%%%%%%%%%%%%%%%%%%% compute weights %%%%%%%%%%%%%%%%%%%%%%
     
    y_k_hat=H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]);    
    theta_hat=bearing_generate(y_k_hat(1),y_k_hat(2),0);% bearing of  y_k_hat

    S_k=H(:,:,i)*P_k_k_1*H(:,:,i)'+Qtm;
    
    if isempty(find(eig(S_k)<=0))==0  % 非正定
        lossnum=lossnum+1;
        lossnum_SRF=lossnum_SRF+1;
        X_k=X_k_k_1;
        P_k=P_k_k_1;
        lossflag=1;
        break;
    end
%     n=mk;
%     bk=bk0;
    if i<=3
        bk=bk0;
        n=mk;
    end
    
    % test if measurements falling in the gate,and extract the measurements
    % in the gate and put it into vector bk
    if i>3
        bk=bk0;
        n=mk;
    end
        
    for j=1:n
    
    
        if i>3
            S_k_sqrt=chol(S_k,'lower'); % S_k_sqrt'*S_k_sqrt=S_k
            S_k_sqrt_inv=inv(S_k_sqrt);
            S_k_sqrt_inv=chol(inv(S_k),'lower');
            s11= S_k_sqrt_inv(1,1);
            s12= S_k_sqrt_inv(1,2);
            s21= S_k_sqrt_inv(2,1);
            s22= S_k_sqrt_inv(2,2);


            g_theta=1/sqrt(bk(:,j)'/S_k*bk(:,j));
            h_theta=bearing_generate(s11*bk(1,j)+s12*bk(2,j),s21*bk(1,j)+s22*bk(2,j),0);
            z=y_k_hat'*S_k_sqrt_inv*[sin(h_theta) cos(h_theta)]';
            f_theta=g_theta*g_theta/sqrt(det(S_k))*1/2/pi*exp(-0.5*(y_k_hat'/S_k*y_k_hat-z*z))*(exp(-0.5*z*z)+sqrt(2*pi)*z*0.5*erfc(-z/sqrt(2)));

            c=1/((1-pc)*f_theta+pc/Vol); %normalizing factor
            w0=c*(1-pc)*f_theta;
            w1=c*pc/Vol;
        else
            w0=1;
            w1=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % %%%%%%%%%%%%%%%%%%%%  compute X_k_0 %%%%%%%%%%%%%%%%%
        %correction steps
        K_k=P_k_k_1*H(:,:,i)'/S_k;
        epsilon_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*bk(:,j)'/S_k*(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]));

        %compute rou2(epsilon_k)
        Fnormal=0.5*erfc(-1*epsilon_k/sqrt(2));
        rou2=(epsilon_k*exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k*epsilon_k+1)*Fnormal)/(exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k)*Fnormal);

        gamma_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*rou2;
        delta_k=1/(bk(:,j)'/S_k*bk(:,j))*(2+epsilon_k*rou2-rou2*rou2);

        %estimation
        X_k_0=X_k_k_1-K_k*H(:,:,i)*X_k_k_1-K_k*utm([2*i-1,2*i])+gamma_k*K_k*bk(:,j);
        P_k_0=P_k_k_1-K_k*H(:,:,i)*P_k_k_1+delta_k*K_k*bk(:,j)*bk(:,j)'*K_k';
                                                                                                                                                                                           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % fusion

        X_k=w0*X_k_0+w1*X_k_k_1;
        P_k=w0*(P_k_0+(X_k-X_k_0)*(X_k-X_k_0)')+w1*(P_k_k_1+(X_k-X_k_k_1)*(X_k-X_k_k_1)');

        X_k_k_1=X_k;
        P_k_k_1=P_k;
    end
    lossflag=0;
end

    