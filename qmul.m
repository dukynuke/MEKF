function [ q_r ] = qmul(q_f, q_s)
%Function used for quaternion multiplication
q_r=[q_s(4) -q_s(3) q_s(2) q_s(1)
    q_s(3) q_s(4) -q_s(1) q_s(2)
    -q_s(2) q_s(1) q_s(4) q_s(3)
    -q_s(1) -q_s(2) -q_s(3) q_s(4)]*q_f;

end

