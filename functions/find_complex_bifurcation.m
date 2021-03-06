function f = find_complex_bifurcation(x, p1, p3, p4, p5, p6)
    if p3 == 0
        exp_x2 = exp(x(2));
    else
        exp_x2 = exp(x(2) / (1 + x(2) / p3));
    end
    f = [-p1 * x(1) + x(3) * (1 - x(1)) * exp_x2;...
         -p1 * x(2) + x(3) * p4 * (1 - x(1)) * exp_x2 - p5 * (x(2) - p6);...
         df1x1(x(1), x(2), p1, x(3), p3, p4, p5, p6)... 
         + df2x2(x(1), x(2), p1, x(3), p3, p4, p5, p6)];
end