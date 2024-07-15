%%matlab R2024a 
%@author: WANG YIMING
%@mail: wang.yiming.36j@st.kyoto-u.ac.jp
%@date: 2024/07/05


x_step = 0.01;
x = 0:x_step:4;
t = [0 0 0 0 4 4 4 4];
degrees = 3;

y = zeros(length(x), length(t) - (degrees + 1));
start_flag = 1;
end_flag = 1;

%to match the length of X and B
%set start point
for i=start_flag:length(x)
    if x(i) >= t(1)
        start_flag = i;
        break;
    end
end

%set end point
for i=start_flag:length(x)
    if x(i) >= t(length(t))
        end_flag = i;
        break;
    end
end

y(start_flag:end_flag, :) = BSplineGenerate(t, degrees, x_step);

figure;
hold on;
for i=1:size(y,2)
    plot(x, y(:,i));
end
hold off;
xlabel('X');
ylabel('B');
title('BSpline');

%B-Spline count
%x is the variable
%degrees is B's degree
%t is note column
function result=BSplineGenerate(t, degrees, x_step)
    t_x = t(1):x_step:t(length(t));
    B = zeros(length(t_x), length(t) - (degrees + 1));
    for i=1:(size(t, 2) - (degrees + 1))
        for k=0:degrees
            for index=1:length(t_x)
                if t_x(index)==t(i+k)
                    countStart = index;
                end
                if t_x(index)==t(i+k+1)
                    countend = index;
                    break;
                end
            end
            pmt = BSplineCoefGenerate(k, degrees);
            %B_i generate
            B(countStart:countend, i) = BSplineCount(i, k ,t, degrees, t_x(countStart:countend), pmt);
        end
    end
    result = B;
end

%the coef. of b-spline generate.
%like 000, 001, 010, 100, ...
%k-- as 0 effcting
%i++, k-- as 1 effcting
function result=BSplineCoefGenerate(k, degrees)
    array=zeros(1,degrees);
    if k>0
        for j=1:k
            array(j)=1;
        end
    end
    result = unique(perms(array), "rows");
end

%B-Spline count
%return [:, 1]
%b(_i^k)=(x - t_i)/(t_i+1 - t_i) => 
%1 / (t_i+1 - t_i) * (x - t_i)
function result=BSplineCount(i, k, t, degrees, t_x, coefNumber)
    B_ipkToNext = zeros(size(t_x, 2), 1);

    [coef, b] = BSplineCoefCount(i, t, degrees, coefNumber);

    for j=1:length(B_ipkToNext)
        %uniform: x<t(i+k+1)
        %open uniform: x<=t(i+k+1)
        if t(i+k)<=t_x(j) && t_x(j)<=t(i+k+1)
            for m=1:size(coefNumber, 1)
                Bn = coef(m,1)*(t_x(j)-b(m,1));
                for n=2:degrees
                    Bn = Bn*coef(m,n)*(t_x(j)-b(m,n));
                end
                B_ipkToNext(j) = B_ipkToNext(j) + Bn;
            end
        else
            B_ipkToNext(j)=0;
        end
    end

    result = B_ipkToNext;
end

function [coef, b]=BSplineCoefCount(i, t, degrees, coefNumber)
    coef = zeros(size(coefNumber));
    b = zeros(size(coefNumber));
    for m=1:size(coefNumber, 1)
        idx = i;
        kdx = degrees;
        for n=1:degrees
            if coefNumber(m,n) == 0
                if (t(idx+1) - t(idx))~=0
                    coef(m, n) = 1 / (t(idx+1) - t(idx));
                    b(m, n) = t(idx);
                end
                kdx = kdx - 1;
            end
            if coefNumber(m,n) ==1
                if (t(idx+kdx+1) - t(idx+1))~=0
                    coef(m, n) = -1 / (t(idx+kdx+1) - t(idx+1));
                    b(m, n) = t(idx+kdx+1);
                end
                idx = idx + 1;
                kdx = kdx - 1;
            end
        end
    end
end

