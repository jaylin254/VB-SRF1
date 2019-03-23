%由于atan函数的值在-pi/2 到pi/2之间，不符合实际的角度值，再求该角度的sin和cos函数时
%会带来计算符号的错误，所以在此进行处理

function angle=bearing_generate(x,y,sigma)

% sigma：standard deviation of the bearing
    if x>=0&&y>0
            angle=atan(x/y)+randn*sigma;   % measurementsatan
    end
    if y<0&&x>=0
       angle=pi+atan(x/y)+randn*sigma;   % measurementsatan

    end
    if y<0&&x<=0
        angle=-pi+atan(x/y)+randn*sigma;   % measurementsatan

    end
    if y>0&&x<=0
        angle=atan(x/y)+randn*sigma;   % measurementsatan

    end
end