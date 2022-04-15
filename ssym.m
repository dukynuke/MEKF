function [ mat_out ] = ssym(vec_in)
%Function used to build skew-symmetric matrix
mat_out=[0 -vec_in(3) vec_in(2)
    vec_in(3) 0 -vec_in(1)
    -vec_in(2) vec_in(1) 0];

end