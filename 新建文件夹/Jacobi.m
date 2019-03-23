function [H,hh]=Jacobi(X_k_k_1,utm,i)
% compute H for angle-only measurement
% i denotes the ith angle measurement from the target to the sonobuoys    
H=zeros(1,12);
if i==1
    rh_squa=(utm(1)-X_k_k_1(5))^2+(utm(2)-X_k_k_1(6))^2;

    H(1,5)=-(utm(2)-X_k_k_1(6))/rh_squa;
    H(1,6)=(utm(1)-X_k_k_1(5))/rh_squa; 
    
    hh(1)=bearing_generate(utm(1)-X_k_k_1(5),utm(2)-X_k_k_1(6),0); 
   
end
if i==2
    rh_squa=(utm(1)-X_k_k_1(7))^2+(utm(2)-X_k_k_1(8))^2;

    H(1,7)=-(utm(2)-X_k_k_1(7))/rh_squa;
    H(1,8)=(utm(1)-X_k_k_1(8))/rh_squa; 
    
    hh(1)=bearing_generate(utm(1)-X_k_k_1(7),utm(2)-X_k_k_1(8),0); 
   
end
if i==3
    rh_squa=(utm(1)-X_k_k_1(9))^2+(utm(2)-X_k_k_1(10))^2;

    H(1,9)=-(utm(2)-X_k_k_1(9))/rh_squa;
    H(1,10)=(utm(1)-X_k_k_1(10))/rh_squa; 
    
    hh(1)=bearing_generate(utm(1)-X_k_k_1(9),utm(2)-X_k_k_1(10),0); 
   
end


if i==4
    rh_squa=(X_k_k_1(1)-X_k_k_1(5))^2+(X_k_k_1(3)-X_k_k_1(6))^2;

    H(1,1)=(X_k_k_1(3)-X_k_k_1(6))/rh_squa;
    H(1,3)=-(X_k_k_1(1)-X_k_k_1(5))/rh_squa;
    H(1,5)=-(X_k_k_1(3)-X_k_k_1(6))/rh_squa;
    H(1,6)=(X_k_k_1(1)-X_k_k_1(5))/rh_squa;
    
    %%% compute hh %%%%%%%%%%%%%%%%%
    hh(1)=bearing_generate((X_k_k_1(1)-X_k_k_1(5)),(X_k_k_1(3)-X_k_k_1(6)),0);

end

if i==5
    rh_squa=(X_k_k_1(1)-X_k_k_1(7))^2+(X_k_k_1(3)-X_k_k_1(8))^2;

    H(1,1)=(X_k_k_1(3)-X_k_k_1(8))/rh_squa;
    H(1,3)=-(X_k_k_1(1)-X_k_k_1(7))/rh_squa;
    H(1,7)=-(X_k_k_1(3)-X_k_k_1(8))/rh_squa;
    H(1,8)=(X_k_k_1(1)-X_k_k_1(7))/rh_squa;
    
    %%% compute hh %%%%%%%%%%%%%%%%%
    hh(1)=bearing_generate((X_k_k_1(1)-X_k_k_1(7)),(X_k_k_1(3)-X_k_k_1(8)),0);

end

if i==6
    rh_squa=(X_k_k_1(1)-X_k_k_1(9))^2+(X_k_k_1(3)-X_k_k_1(10))^2;

    H(1,1)=(X_k_k_1(3)-X_k_k_1(10))/rh_squa;
    H(1,3)=-(X_k_k_1(1)-X_k_k_1(9))/rh_squa;
    H(1,9)=-(X_k_k_1(3)-X_k_k_1(10))/rh_squa;
    H(1,10)=(X_k_k_1(1)-X_k_k_1(9))/rh_squa;
    
    %%% compute hh %%%%%%%%%%%%%%%%%
    hh(1)=bearing_generate((X_k_k_1(1)-X_k_k_1(9)),(X_k_k_1(3)-X_k_k_1(10)),0);

end