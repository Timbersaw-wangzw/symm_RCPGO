function J = Jl(vec)
% leftJacobian of SO3 and SE3. vec is 3x1 or 6x1 
len=length(vec);
switch len
    case 3
        angle=norm(vec);
        if(angle > 1e-6)
            B = (1 - cos(angle))/(angle*angle);
            C = (angle - sin(angle))/(angle*angle*angle);
        else
            angle2 = angle * angle;
            B = 0.5 - angle2/24 + angle2 * angle2/720;
            C = 1/6 - angle2/120;
        end
        vHat = skew(vec);
        J= eye(3) + B * vHat + C * vHat * vHat ;
    case 6
        J=zeros(6,6);
        t_Jlw= Jl(vec(4:6));
        J(1:3,1:3) = t_Jlw;
        J(1:3,4:6)=Ql(vec);
        J(4:6,4:6)=t_Jlw;
end
end

