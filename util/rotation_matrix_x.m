function R = rotation_matrix_x(th)
R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];