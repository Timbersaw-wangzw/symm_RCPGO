function M=InvJl(v)
len=length(v);
M=zeros(3,3);
switch len
    case 3
        angle=norm(v);
        if(angle > 1e-6)
            D = (1/angle - sin(angle)/(2*(1-cos(angle))))/angle;
        
        else
            angle2 = angle * angle;
            D = ( 1/12 - angle2/180 ) / ( 1 - angle2/12 + (angle2 * angle2)/360 );
        end
        vHat = skew(v);
        M= eye(3) - 0.5 * vHat + D * vHat * vHat;
    case 6
        M=zeros(6,6);
        t_InvJlw= InvJl(v(4:6));
        M(1:3,1:3)=t_InvJlw;
        M(1:3,4:6)=-t_InvJlw*Ql(v)*t_InvJlw;
        M(4:6,4:6)=t_InvJlw;
end
end