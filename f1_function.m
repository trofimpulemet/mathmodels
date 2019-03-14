function f = f1_function(x1, x2, p2, p1, p3, p4, p5, p6)
    if p3 == 0
        exp_x2 = exp(x2);
    else
        exp_x2 = exp(x2 / (1 + x2 / p3));
    end
    f = -p1 * x1 + p2 * (1 - x1) * exp_x2;
end

