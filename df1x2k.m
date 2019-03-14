function f = df1x2k(x)
    p1 = 1; 
    p3 = 20;
    p4 = 10;
    p5 = 0.6;
    p6 = -5;
%     if p3 == 0
%         exp_x2 = exp(x2);
%         exp_arg = 1;
%     else
    exp_arg = 1 + x / p3;
    x1 = (x * (p1 + p5) - p5 * p6) / (p1 * p4);
    p2_exp_x2 = (p1 * x1) / (1 - x1);
    f = (p2_exp_x2 * (1 - x1)) / (exp_arg * exp_arg);
end