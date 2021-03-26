clear
clc

% Points in data
iteration_time = 240*3;
iteration_total = 1;
total_time = iteration_time*iteration_total;

% grid dimensions
x_dim = 10;
y_dim = 10;
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

x_it_temp = x_it;
x_it = [x_it y_it];
y_it = [y_it x_it_temp];

filename = 'shelves_rotation';

save(filename,'x_unique','y_unique','x_it','y_it')