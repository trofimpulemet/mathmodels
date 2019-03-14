function f = df1x1(x1, x2, p1, p2, p3, p4, p5, p6)
    if p3 == 0
        exp_x2 = exp(x2);
    else
        exp_x2 = exp(x2 / (1 + x2 / p3));
    end
    p2_exp_x2 = (x2 * (p1 + p5) - p5 * p6) / (p1 - p4 - x2 * (p1 + p5) - p5  * p6);
    f = -p1 - p2 * exp_x2; 
end
