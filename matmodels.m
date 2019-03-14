clc;
clear;

% p1 = 1; %а
% p3 = 20;
% p4 = 10;
% p5 = 0.6;
% p6 = -5;

p1 = 0.5; %а
p3 = 0;
p4 = 14; % 8, 10, 12, 14
p5 = 0.8;
p6 = 0;

eps = 1e-10;

x2_begin = -2;
x2_end = 5;
h = 0.1;
x2_interval = zeros(1, ((x2_end - x2_begin) / h) + 1);
x2_interval(1) = x2_begin;
for i = 1:length(x2_interval)
    x2_interval(i + 1) = x2_interval(i) + h;
end
m = 1;
for i = 1:length(x2_interval)
    x1_check = x1_function(x2_interval(i), p1, p4, p5, p6);
    p2_check = p2_function(x2_interval(i), p1, p3, p4, p5, p6);
    if x1_check > eps && x1_check < 1 && p2_check > 0
        x1(m) = x1_check;
        p2(m) = p2_check;
        x2(m) = x2_interval(i);
        m = m + 1;
    end
end

lambda1 = zeros(1, length(x1));
lambda2 = zeros(1, length(x1));
for i = 1: length(x1)
    a11 = df1x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a12 = df1x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a21 = df2x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a22 = df2x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
 
    A = [a11, a12; a21, a22]; 
    eig_A = eig(A);
    lambda1(i) = eig_A(1);
    lambda2(i) = eig_A(2);
end

C(length(x1), 5) = zeros;
for i = 1:length(x1)
    C(i, 1) = x1(i);
    C(i, 2) = x2(i);
    C(i, 3) = p2(i);
    C(i, 4) = lambda1(i);
    C(i, 5) = lambda2(i);  
end
T = array2table(C, 'VariableNames', {'x1', 'x2', 'p2', 'l1', 'l2'});

m = 1;
for i = 1: length(x1) - 1 
    a11 = df1x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a12 = df1x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a21 = df2x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a22 = df2x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    
    b11 = df1x1(x1(i + 1), x2(i + 1), p1, p2(i + 1), p3, p4, p5, p6);
    b12 = df1x2(x1(i + 1), x2(i + 1), p1, p2(i + 1), p3, p4, p5, p6);
    b21 = df2x1(x1(i + 1), x2(i + 1), p1, p2(i + 1), p3, p4, p5, p6);
    b22 = df2x2(x1(i + 1), x2(i + 1), p1, p2(i + 1), p3, p4, p5, p6);
    
    A = [a11, a12; a21, a22]; 
    ei1 = eig(A); 
    
    B = [b11, b12; b21, b22]; 
    ei2 = eig(B);

    if ((real(ei1(1)) < 0 && real(ei1(2)) < 0)... 
       && (real(ei2(1)) >= 0 || real(ei2(2)) >= 0))...
       || ((real(ei1(1)) >= 0 || real(ei1(2)) >= 0)...
       && (real(ei2(1)) < 0 && real(ei2(2)) < 0))   
        x1trans(m) = x1(i);
        x2trans(m) = x2(i);
        p2trans(m) = p2(i);
%         v(m) = i;
        m = m + 1;
    end
end

check_exist_trans_point = exist('x1trans', 'var');
if check_exist_trans_point == 1
    index_real = 1;
    index_complex = 1;
    for i = 1:length(x1trans)
        options = optimset('TolFun',1.0e-12);
        x0 = [x1trans(i), x2trans(i), p2trans(i)]; 
        xr_complex = fsolve(@find_complex_bifurcation, x0, options, p1, p3, p4, p5, p6);
        xr_real = fsolve(@find_real_bifurcation, x0, options, p1, p3, p4, p5, p6);

        a11 = df1x1(xr_complex(1), xr_complex(2), p1, xr_complex(3), p3, p4, p5, p6);
        a12 = df1x2(xr_complex(1), xr_complex(2), p1, xr_complex(3), p3, p4, p5, p6);
        a21 = df2x1(xr_complex(1), xr_complex(2), p1, xr_complex(3), p3, p4, p5, p6);
        a22 = df2x2(xr_complex(1), xr_complex(2), p1, xr_complex(3), p3, p4, p5, p6);

        b11 = df1x1(xr_real(1), xr_real(2), p1, xr_real(3), p3, p4, p5, p6);
        b12 = df1x2(xr_real(1), xr_real(2), p1, xr_real(3), p3, p4, p5, p6);
        b21 = df2x1(xr_real(1), xr_real(2), p1, xr_real(3), p3, p4, p5, p6);
        b22 = df2x2(xr_real(1), xr_real(2), p1, xr_real(3), p3, p4, p5, p6);

        A = [a11, a12; a21, a22];
        eig_complex = eig(A);

        B = [b11, b12; b21, b22];
        eig_real = eig(B);

         if (real(eig_complex(1)) < eps && real(eig_complex(2)) < eps)...
            && (abs(eig_complex(1)) == abs(eig_complex(2)))...
            && (imag(eig_complex(1)) ~= 0 && imag(eig_complex(2)) ~= 0)
            x1_complex_bf(index_complex) = xr_complex(1);
            x2_complex_bf(index_complex) = xr_complex(2);
            p2_complex_bf(index_complex) = xr_complex(3);
            index_complex = index_complex + 1;
         end

         if eig_real(1) < eps || eig_real(2) < eps
             x1_real_bf(index_real) = xr_real(1);
             x2_real_bf(index_real) = xr_real(2);
             p2_real_bf(index_real) = xr_real(3);
         index_real = index_real + 1;
         end
    end
end

% Построение графиков

% x1(p2)
f1 = figure('Name', 'x1(p2)');
grid on;
title('x1(p2)');
xlabel('p2');
ylabel('x1');
hold on;
plot(p2, x1, 'Marker', 'none', 'Color', [0, 0, 0], 'LineWidth', 1.5);
set(gca, 'FontSize', 24)

for i = 1: length(x1)
    
    a11 = df1x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a12 = df1x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a21 = df2x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a22 = df2x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
 
    A = [a11, a12; a21, a22]; 
    eig_A = eig(A); 
   if (real(eig_A) < 0) 
       h1 = plot(p2(i), x1(i), 'Marker', '.', 'MarkerSize', 14, 'Color', [0, 0.8, 0]);
   else
       h2 = plot(p2(i), x1(i), 'Marker', '.', 'MarkerSize', 14, 'Color', 'r');
   end
end


% Наложение на график точек бифуркации
check_existing_real_bf = exist('x1_real_bf', 'var');
if check_existing_real_bf == 1
    for i = 1:length(x1_real_bf)
        h3 = plot(p2_real_bf(i), x1_real_bf(i), 'Marker', '.', 'MarkerSize', 32, 'Color', [0, 0, 0]);
    end
end

check_existing_complex_bf = exist('x1_complex_bf', 'var');
if check_existing_complex_bf == 1
    for i = 1:length(x1_complex_bf)
        h4 = plot(p2_complex_bf(i), x1_complex_bf(i), 'Marker', '.', 'MarkerSize', 32, 'Color', [0, 0, 1]);
    end
end

if check_existing_real_bf == 1 && check_existing_complex_bf == 1
    l1 = legend([h1, h2, h3, h4], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                 'Точка вещественной бифуркации', 'Точка комплексной бифуркации');
    l1 = legend('boxoff');
    l1 = legend('Location', 'best');
    l1.FontSize = 22;
elseif check_existing_real_bf == 1 && check_existing_complex_bf == 0
    l1 = legend([h1, h2, h3], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                'Точка вещественной бифуркации');
    l1 = legend('boxoff');
    l1 = legend('Location', 'best');
    l1.FontSize = 22;
elseif check_existing_real_bf == 0 && check_existing_complex_bf == 1
    l1 = legend([h1, h2, h4], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                'Точка комплексной бифуркации');
    l1 = legend('boxoff');
    l1 = legend('Location', 'best');
    l1.FontSize = 22;
elseif check_existing_real_bf == 0 && check_existing_complex_bf == 0
    check_exist_stable_point = exist('h1', 'var');
    check_exist_unstable_point = exist('h2', 'var');
    if check_exist_stable_point == 1 && check_exist_unstable_point == 1
        l1 = legend([h1, h2], 'Устойчивая точка решения', 'Неустойчивая точка решения');
        l1 = legend('boxoff');
        l1 = legend('Location', 'best');
        l1.FontSize = 22;
    elseif check_exist_stable_point == 1 && check_exist_unstable_point == 0
        l1 = legend(h1, 'Устойчивая точка решения');
        l1 = legend('boxoff');
        l1 = legend('Location', 'best');
        l1.FontSize = 22;
    elseif check_exist_stable_point == 0 && check_exist_unstable_point == 1
        l1 = legend(h2, 'Неустойчивая точка решения');
        l1 = legend('boxoff');
        l1 = legend('Location', 'best');
        l1.FontSize = 22;
    end
end
hold off;

% x2(p2)
f2 = figure('Name', 'x2(p2)');
grid on;
title('x2(p2)');
xlabel('p2');
ylabel('x2');
hold all;
plot(p2, x2, 'Marker', 'none', 'Color', [0, 0, 0], 'LineWidth', 1.5);
set(gca, 'FontSize', 24)

for i = 1: length(x1)
    
    a11 = df1x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a12 = df1x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a21 = df2x1(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
    a22 = df2x2(x1(i), x2(i), p1, p2(i), p3, p4, p5, p6);
 
    A = [a11, a12; a21, a22]; 
    eig_A = eig(A); 
   if (real(eig_A) < 0) 
       h5 = plot(p2(i), x2(i), 'Marker', '.', 'MarkerSize', 14, 'Color', [0, 0.8, 0]);
   else
       h6 = plot(p2(i), x2(i), 'Marker', '.', 'MarkerSize', 14, 'Color', 'r');
   end
end


% Наложение на график точек бифуркации
check_existing_real_bf = exist('x2_real_bf', 'var');
if check_existing_real_bf == 1
    for i = 1:length(x2_real_bf)
        h7 = plot(p2_real_bf(i), x2_real_bf(i), 'Marker', '.', 'MarkerSize', 32, 'Color', [0, 0, 0]);
    end
end

check_existing_complex_bf = exist('x2_complex_bf', 'var');
if check_existing_complex_bf == 1
    for i = 1:length(x2_complex_bf)
        h8 = plot(p2_complex_bf(i), x2_complex_bf(i), 'Marker', '.', 'MarkerSize', 32, 'Color', [0, 0, 1]);
    end
end

if check_existing_real_bf == 1 && check_existing_complex_bf == 1
    l2 = legend([h5, h6, h7, h8], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                 'Точка вещественной бифуркации', 'Точка комплексной бифуркации');
    l2 = legend('boxoff');
    l2 = legend('Location', 'best');
    l2.FontSize = 22;
elseif check_existing_real_bf == 1 && check_existing_complex_bf == 0
    l2 = legend([h5, h6, h7], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                'Точка вещественной бифуркации');
    l2 = legend('boxoff');
    l2 = legend('Location', 'best');
    l2.FontSize = 22;
elseif check_existing_real_bf == 0 && check_existing_complex_bf == 1
    l2 = legend([h5, h6, h8], 'Устойчивая точка решения', 'Неустойчивая точка решения',...
                'Точка комплексной бифуркации');
    l2 = legend('boxoff');
    l2 = legend('Location', 'best');
    l2.FontSize = 22;
elseif check_existing_real_bf == 0 && check_existing_complex_bf == 0
    check_exist_stable_point = exist('h5', 'var');
    check_exist_unstable_point = exist('h6', 'var');
    if check_exist_stable_point == 1 && check_exist_unstable_point == 1
        l2 = legend([h5, h6], 'Устойчивая точка решения', 'Неустойчивая точка решения');
        l2 = legend('boxoff');
        l2 = legend('Location', 'best');
        l2.FontSize = 22;
    elseif check_exist_stable_point == 1 && check_exist_unstable_point == 0
        l2 = legend(h5, 'Устойчивая точка решения');
        l2 = legend('boxoff');
        l2 = legend('Location', 'best');
        l2.FontSize = 22;
    elseif check_exist_stable_point == 0 && check_exist_unstable_point == 1
        l2 = legend(h6, 'Неустойчивая точка решения');
        l2 = legend('boxoff');
        l2 = legend('Location', 'best');
        l2.FontSize = 22;
    end
end
hold off;

% Блок проверки
check_right = zeros(length(x1), 5);
for i = 1:length(x1)
    check_right(i, 1) = x1(i);
    check_right(i, 2) = x2(i);
    check_right(i, 3) = p2(i);
    check_right(i, 4) = f1_function(x1(i), x2(i), p2(i), p1, p3, p4, p5, p6);
    check_right(i, 5) = f2_function(x1(i), x2(i), p2(i), p1, p3, p4, p5, p6);
end
T_check = array2table(check_right);
xlswrite('check.xlsx', check_right, 'A2:E63');

    































