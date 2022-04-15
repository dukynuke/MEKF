
function dcm= quat_2_dcm(quat)

q1 = quat(1);
q2 = quat(2);
q3 = quat(3);
q4 = quat(4);


dcm(1,1)=q1^2-q2^2-q3^2+q4^2;
dcm(1,2)=2*(q1*q2+q3*q4);
dcm(1,3)=2*(q1*q3-q2*q4);

dcm(2,1)=2*(q1*q2-q3*q4);
dcm(2,2)=(q2^2-q1^2-q3^2+q4^2);
dcm(2,3)=2*(q2*q3+q1*q4);

dcm(3,1)=2*(q1*q3+q2*q4);
dcm(3,2)=2*(q2*q3-q1*q4);
dcm(3,3)=(q3^2-q1^2-q2^2+q4^2);