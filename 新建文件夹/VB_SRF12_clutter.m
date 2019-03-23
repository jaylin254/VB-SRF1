function [X_k,P_k,eta_k,alpha1_k,alpha2_k,lossflag]=VB_SRF12_clutter(X0,P0,Z_k,utm,eta_k_1,alpha1_k_1,alpha2_k_1)
% 利用VB方法来估计杂波出现概率

global h q1 q2 q3 sigma sigma_prime  lossnum lossnum_VB
global B
% shifted rayleigh filter
% model : X_k=Fai*X_k_1+V_k
% X=12 dimentional
% z_k1 is the conditional angle ,like bearing
% R_k is the variance of z_k1
rou=1-exp(-4);

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
M=2;
% projected measurement 
% bk=[sin(Z_k(1)),cos(Z_k(1)),sin(Z_k(2)),cos(Z_k(2)),sin(Z_k(3)),cos(Z_k(3)),sin(Z_k(4)),cos(Z_k(4)),sin(Z_k(5)),cos(Z_k(5)),sin(Z_k(6)),cos(Z_k(6))]';


X_k_1=X0;
P_k_1 = P0; % initial estimation covariance

%filter

%prediction
X_k_k_1=Fai*X_k_1;
P_k_k_1=Fai*P_k_1*Fai'+Qs;
eta_k_k_1=rou*eta_k_1;
alpha2_k_k_1=rou*alpha2_k_1;
alpha1_k_k_1=rou*alpha1_k_1;

%update

for i=1:6
        
    bk0=[sin(Z_k{i}),cos(Z_k{i})]';
    mk=size(bk0,2); % number of measurements for one sensor
        
    if i<=3
        Qtm=Qtr+(norm(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]),2).^2+trace(H(:,:,i)*P_k_k_1*H(:,:,i)'))*sigma*sigma*eye(2);
    else
        Qtm=Qtr+(norm(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]),2).^2+trace(H(:,:,i)*P_k_k_1*H(:,:,i)'))*sigma_prime*sigma_prime*eye(2);
    end
    
    S_k=H(:,:,i)*P_k_k_1*H(:,:,i)'+Qtm;
   
    Qtm_sqrt_inv=chol(inv(Qtm),'lower');
    qt11= Qtm_sqrt_inv(1,1);
    qt12= Qtm_sqrt_inv(1,2);
    qt21= Qtm_sqrt_inv(2,1);
    qt22= Qtm_sqrt_inv(2,2);
   
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
        
        if i<=3
            %correction steps
            K_k=P_k_k_1*H(:,:,i)'/S_k;
            epsilon_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*bk(:,j)'/S_k*(H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]));

            %compute rou2(epsilon_k)
            Fnormal=0.5*erfc(-1*epsilon_k/sqrt(2));
            rou2=(epsilon_k*exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k*epsilon_k+1)*Fnormal)/(exp(-0.5*epsilon_k*epsilon_k)+sqrt(2*pi)*(epsilon_k)*Fnormal);

            gamma_k=(bk(:,j)'/S_k*bk(:,j))^(-0.5)*rou2;
            delta_k=1/(bk(:,j)'/S_k*bk(:,j))*(2+epsilon_k*rou2-rou2*rou2);

            %estimation
            X_k=X_k_k_1-K_k*H(:,:,i)*X_k_k_1-K_k*utm([2*i-1,2*i])+gamma_k*K_k*bk(:,j);
            P_k=P_k_k_1-K_k*H(:,:,i)*P_k_k_1+delta_k*K_k*bk(:,j)*bk(:,j)'*K_k';
            
            X_k_k_1=X_k;
            P_k_k_1=P_k;
            
        end
        
        if i>3
            
            X_k=X_k_k_1;
            P_k=P_k_k_1;
            eta_k(i-3)=eta_k_k_1(i-3);
            alpha1_k(i-3)=alpha1_k_k_1(i-3);
            alpha2_k(i-3)=alpha2_k_k_1(i-3);
                          
            % %%%%%%%%%%%%%%%%%%%%%%% compute the density function of theta %%%%%%%%%%%%%%%%%%%%%%
     
%             S_k_sqrt=chol(S_k,'lower'); % S_k_sqrt'*S_k_sqrt=S_k
%             S_k_sqrt_inv=inv(S_k_sqrt);
            if isempty(find(eig(S_k)<=0))==0  % 非正定
                    lossnum=lossnum+1;
                    lossnum_VB=lossnum_VB+1;
                    X_k=X_k_k_1;
                    P_k=P_k_k_1;
                    lossflag=1;
                    break;
            end
            S_k_sqrt_inv=chol(inv(S_k),'lower');
            s11= S_k_sqrt_inv(1,1);
            s12= S_k_sqrt_inv(1,2);
            s21= S_k_sqrt_inv(2,1);
            s22= S_k_sqrt_inv(2,2);
            
            y_k_hat=H(:,:,i)*X_k_k_1+utm([2*i-1,2*i]);  
            g_theta=1/sqrt(bk(:,j)'/S_k*bk(:,j));
            h_theta=bearing_generate(s11*bk(1,j)+s12*bk(2,j),s21*bk(1,j)+s22*bk(2,j),0);
            z=y_k_hat'*S_k_sqrt_inv*[sin(h_theta) cos(h_theta)]';
            f_theta=g_theta*g_theta/sqrt(det(S_k))*1/2/pi*exp(-0.5*(y_k_hat'/S_k*y_k_hat-z*z))*(exp(-0.5*z*z)+sqrt(2*pi)*z*0.5*erfc(-z/sqrt(2)));
           
            for m=1:M
                
               %compute weights
                c=1/((1-eta_k(i-3))*f_theta+eta_k(i-3)/Vol); %normalizing factor
                w0=c*(1-eta_k(i-3))*f_theta;
                w1=c*eta_k(i-3)/Vol;
        
                % %%%%%%%%%%%%%%%%%%%%  compute X_k %%%%%%%%%%%%%%%%%
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
             
                % fusion
                X_k=w0*X_k_0+w1*X_k_k_1;
                P_k=w0*(P_k_0+(X_k-X_k_0)*(X_k-X_k_0)')+w1*(P_k_k_1+(X_k-X_k_k_1)*(X_k-X_k_k_1)');
                
                % %%%%%%%%%%%%%% update eta_k %%%%%%%%%%%%%%%%%%%%%%%%%
               
                % first compute ln f(theta|xk)
                
                
                y_k_hat_xn=H(:,:,i)*X_k+utm([2*i-1,2*i]);   
                g_theta=1/sqrt(bk(:,j)'/Qtm*bk(:,j));
                h_theta=bearing_generate(qt11*bk(1,j)+qt12*bk(2,j),qt21*bk(1,j)+qt22*bk(2,j),0);
                z=y_k_hat_xn'*Qtm_sqrt_inv*[sin(h_theta) cos(h_theta)]';
                f_theta_xn=g_theta*g_theta/sqrt(det(Qtm))*1/2/pi*exp(-0.5*(y_k_hat_xn'/Qtm*y_k_hat_xn-z*z))*(exp(-0.5*z*z)+sqrt(2*pi)*z*0.5*erfc(-z/sqrt(2)));
                
                % compute eta_k
                feva=psi(alpha1_k(i-3))-psi(alpha2_k(i-3))-log(f_theta_xn)+log(1/Vol);
%                 feva=log( eta_k(i-3))-log( 1-eta_k(i-3))-log(f_theta_xn)+log(1/Vol);
                eta_k(i-3)=exp(feva)/(1+exp(feva));
%                
                %%%%%%%%%%%% update alpha %%%%%%%%%%%%
                alpha2_k(i-3)=alpha2_k(i-3) + 1-eta_k(i-3);
                alpha1_k(i-3)=alpha1_k(i-3)+ eta_k(i-3);
            end
            
                X_k_k_1=X_k;
                P_k_k_1=P_k;

        end
    end
end
eta_k=eta_k';
alpha1_k=alpha1_k';
alpha2_k=alpha2_k';
lossflag=0;
