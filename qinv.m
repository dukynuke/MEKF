function [ q_in ] = qinv(q_duz)
%Function used to take inverse of quaternion
q_in=[-q_duz(1)/(q_duz(1)^2+q_duz(2)^2+q_duz(3)^2+q_duz(4)^2)
    -q_duz(2)/(q_duz(1)^2+q_duz(2)^2+q_duz(3)^2+q_duz(4)^2)
    -q_duz(3)/(q_duz(1)^2+q_duz(2)^2+q_duz(3)^2+q_duz(4)^2)
     q_duz(4)/(q_duz(1)^2+q_duz(2)^2+q_duz(3)^2+q_duz(4)^2)];

end

