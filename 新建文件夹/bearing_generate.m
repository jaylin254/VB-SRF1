%����atan������ֵ��-pi/2 ��pi/2֮�䣬������ʵ�ʵĽǶ�ֵ������ýǶȵ�sin��cos����ʱ
%�����������ŵĴ��������ڴ˽��д���

function angle=bearing_generate(x,y,sigma)

% sigma��standard deviation of the bearing
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