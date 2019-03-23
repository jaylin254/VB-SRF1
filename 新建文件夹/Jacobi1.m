function [H,hh]=Jacobi1(X_k_k_1,utm)
% compute H for angle-only measurement
%  
    H=zeros(1,2);

    rh_squa=(X_k_k_1(1)+utm(1))^2+(utm(2))^2;

    H(1,1)=-utm(2)/rh_squa;
       
    
    %%% compute hh %%%%%%%%%%%%%%%%%
    hh(1)=bearing_generate(utm(2),(X_k_k_1(1)+utm(1)),0);

