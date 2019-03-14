function f = df1x1k(x)
    p1 = 1; 
    p3 = 20;
    p4 = 10;
    p5 = 0.6;
    p6 = -5;
%     if p3 == 0
%         exp_x2 = exp(x2);
%     else
%         exp_x2 = exp(x2 / (1 + x2 / p3));
%     end
    x1 = (x * (p1 + p5) - p5 * p6) / (p1 * p4);
    p2_exp_x2 = (p1 * x1) / (1 - x1);
    f = -p1 - p2_exp_x2; 
end
