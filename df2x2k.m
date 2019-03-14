function f = df2x2k(x)
    p1 = 1; 
    p3 = 20;
    p4 = 10;
    p5 = 0.6;
    p6 = -5;
    
    exp_arg = 1 + x / p3;
    x1 = (x * (p1 + p5) - p5 * p6) / (p1 * p4);
    p2_exp_x2 = (p1 * x1) / (1 - x1);
    f = -p1 - p5 + (p2_exp_x2 * p4 * (1 - x1))... 
              / (exp_arg * exp_arg);
end