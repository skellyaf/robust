function [matrix_X] = crossMatrix (vector_in)


assert(length(vector_in)==3,'Input vector to crossMatrix function is not the correct length');

x = vector_in(1);
y = vector_in(2);
z = vector_in(3);

matrix_X = [0 -z y;
            z 0 -x;
            -y x 0];

end
