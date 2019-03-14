function f = p2_function(x2, p1, p3, p4, p5, p6)
    if p3 == 0
        exp_x2 = exp(x2);
    else
        exp_x2 = exp(x2 / (1 + x2 / p3));
    end
    f = (p1 * x1_function(x2, p1, p4, p5, p6))... 
        / ((1 - x1_function(x2, p1, p4, p5, p6)) * exp_x2);
end