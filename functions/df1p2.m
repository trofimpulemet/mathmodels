function f = df1p2(x1, x2, p1, p2, p3, p4, p5, p6)
    if p3 == 0
        exp_x2 = exp(x2);
    else
        exp_x2 = exp(x2 / (1 + x2 / p3));
    end
    f = (1 - x1) * exp_x2;
end
