function x = cg(f, b, x)
Q_x = Q(reshape(x,128,128));
    r = b - Q_x(:);
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
         Q_p = Q(reshape(p,128,128));
        Ap = Q_p(:);
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end

