clear
clc

% Points in data
iteration_time = 240*6;
iteration_total = 1;
total_time = iteration_time*iteration_total;

% grid dimensions
x_dim = 5;
y_dim = 5;
grid_points_unique = x_dim*y_dim;

% x movement
x_pos = 1:x_dim;
x_neg = x_dim:-1:1;
x = [x_pos x_neg];

% y movement
y = zeros(x_dim);

% y repetition
grid_points_actual = grid_points_unique*((length(x)*length(y))/(x_dim*y_dim));
y_rep = floor(total_time/(grid_points_actual));

% y pre-allocate
y_it = zeros(1,y_rep*grid_points_actual);
y_unique = zeros(1,x_dim*y_dim);

% y coords
for i = 1:y_dim
    if i == 1
        y = ones(1,x_dim)*i;
        y_it(1:length(x)*y_rep) = repmat(y(i,:), 1, y_rep*length(x)/x_dim);
        y_unique = y(i,:);
    else
        temp = ones(1,x_dim)*i;
        y(i,1:x_dim) = temp;
        y_it(:,(i-1)*length(x)*y_rep+1:length(x)*y_rep*i) = repmat(y(i,:), 1, y_rep*length(x)/x_dim);
        y_unique(:,(i-1)*x_dim+1:x_dim*i) = y(i,:);
    end
end

% x repetition
x_rep = y_rep*grid_points_actual/length(x);

% x coords
x_it = repmat(x, 1, x_rep);
x_unique = repmat(x_pos, 1, y_dim);

%x_unique = [ 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5];
%y_unique = [ 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5];

%x_it = [ 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5];
%y_it = [ 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5];

filename = 'shelves';

save(filename,'x_unique','y_unique','x_it','y_it')