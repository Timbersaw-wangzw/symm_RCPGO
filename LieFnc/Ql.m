function M=Ql(vec)
    angle=norm(vec(4:6));
    if angle>1e-6
        sina=sin(angle);
        cosa=cos(angle);
        C = (angle - sina)/(angle*angle*angle);
        E = (1 - 0.5 * angle*angle - cosa)/(angle*angle*angle*angle);
        F = E - (6*C-1)/(2*angle*angle);
    else
        angle2=angle^2;
        C = 1/6 - angle2/120;
        E = -1/24 + angle2/720;
        F = -1/60 + angle2/720;
    end
    uHat=skew(vec(1:3));
    wHat=skew(vec(4:6));
    M = 0.5 * uHat...
        + C * ( wHat*uHat + uHat*wHat + wHat*uHat*wHat )...
        - E * ( wHat*wHat*uHat + uHat*wHat*wHat - 3 * wHat*uHat*wHat)...
        - 0.5 * F * ( wHat*uHat*wHat*wHat + wHat*wHat*uHat*wHat );

end