clear global; clc;
format shortG
format compact
%putting coordinates
coordinates = [0 0; 8 0; 13 0; 8 -3; 8 -6];

%putting the connections
conn = [1 2; 2 3; 2 4; 4 5];
mem=size(conn, 1); %no. of members
nodes = size(coordinates, 1); % number of nodes
ndof = 3*nodes; %degrees of freedom

%array for lengths & angles of the members
lengths = zeros(1, mem);
angles = zeros(1, mem);

%loop for storing angles and lengths of members
for i = 1:mem
    a = conn(i, 1);
    b = conn(i, 2);
    dx = coordinates(a, 1) - coordinates(b, 1);
    dy = coordinates(a, 2) - coordinates(b, 2);
    angles(i)=atan(abs(dy/dx)); 
    lengths(i)=sqrt(dx^2+dy^2);
end

A = 8.0e-3; % m2
I = 75e-6; %m4
E = 200e9; % si units

%allocating the global stiffness matrix, displacement and force vectors
K_global = zeros(ndof, ndof);
D = zeros(ndof, 1);
Q = zeros(ndof, 1);

k_store = cell(mem, 1); %this vector stores the local stiffness matrices of the members
idx_store = cell(mem, 1);%this stores the indices associated with the corresponding 
%nodes according the annoted daigram.
T_matrix_store = cell(mem, 1);% this stores the Transformation matrix of the members

%storing the indices associated with the members
idx_store{1} = [13 14 15 1 2 3];
idx_store{2} = [1 2 3 10 11 12];
idx_store{3} = [1 2 3 4 5 6];
idx_store{4} = [4 5 6 7 8 9];

%calculating local stiffness matrix of members as well as global stiffness
for i = 1:mem
    idx = idx_store{i}; % iterating though each element in turn
    
    L= lengths(i);
    theta = angles(i);
    S = sin(theta);
    C = cos(theta);
    
    T = [C, S, 0, 0, 0, 0;
        -S, C, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, C, S, 0;
        0, 0, 0, -S, C, 0;
        0, 0, 0, 0, 0, 1;];
    
    T_matrix_store{i} = T;
    
    c1 = A*E/L;
    c2 = E*I/(L^3);
    
    %IN CASE OF 2EI WE ASSUMED I AS DOUBLE INSTEAD OF E.
    %member 1 has 2EI others have EI
    if i==1
        c2 = 2*E*I/(L^3);
    end
    
    % WE divided the k(local member stiffness matrix) into 2 parts. Beam and Bar,
    % The bar consists of all the terms of AE/L and beam contains the ,
    % EI/L2 parts and later add both of them.
    
    k_bar = c1*[1, 0, 0, -1, 0, 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;
                -1, 0, 0, 1, 0, 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;];
            
    k_beam = c2*[0,  0,   0, 0,   0,   0;
                 0, 12, 6*L, 0, -12, 6*L;
                 0, 6*L, 4*L^2, 0, -6*L, 2*L^2;
                 0,   0,    0, 0,  0,    0;
                 0, -12, -6*L, 0, 12, -6*L;
                 0, 6*L, 2*L^2, 0, -6*L, 4*L^2;];
             
    k = k_bar + k_beam; %local member stiffness matrix
    k_store{i} = k; %storing the local member matrix in k_store
    
    K = T'*k*T;     % global member stiffness matrix
    %correspondingly adding the indices to form the gloabal Stiffness
    %matirx
    K_global(idx, idx) = K_global(idx, idx) + K;
    
end

%the following boolean is declared to segregate K11 and K21 from GSM
free_dof = logical([1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]); % there are total 15 degrees of freedom,
%out of which first 6 are free and unknown marked as 1 where as rest are
%marked 0

constrained_dof=not(free_dof);%negation of free_dof

%known forces and moments as stated in problem and FEMs
Q(2) = -16.25e3;
Q(3) = -14.58e3;
Q(4) = 20e3;

K11 = K_global(free_dof, free_dof);
K21 = K_global(constrained_dof, free_dof);

%calculating unknown displacements;
D(free_dof) = K11\Q(free_dof); % since the known displacement is zero, the zero vector is not added.

%calculating unknown nodal forces;
Q(constrained_dof) = K21*D(free_dof);

disp("Geometrical and Material Parameters:")
fprintf("Area of cross-section: %f m^2\n", A); 
fprintf("Second Moment of Inertia: %f m^4\n", I); 
fprintf("Young's Modulus of Elasticity: %f GPa\n\n", 200); 

disp("the Displacement/Angle and Force/Moment of the indices are:(respectively)all in SI units");

for i = 1:ndof
   fprintf("%d : %f, %f\n", i, D(i), Q(i));
end

%internal loadings of members
q = cell(mem, 1); %vector for storing the q values of members

for i = 1:mem
    idx = idx_store{i};
    q{i} = k_store{i}*T_matrix_store{i}*D(idx);
end

fprintf("\n\n");

%since the member 2 is having trapeziodal UDL we will superpose the FEM and
% sheer forces with the internal loadings

q{2}(2) = q{2}(2)+16250;
q{2}(5) = q{2}(5)+21250;
q{2}(6) = q{2}(6)-16670;
q{2}(3) = q{2}(3)+14580;

disp("The internal loadings of the members are:(all in SI units)");
for i = 1:mem
    fprintf("Member %d:\n", i);
    for j = 1:6
       fprintf("q %d = %f\n",idx_store{i}(j), q{i}(j));
    end
    fprintf("\n\n");
end

%all the measurements are in SI units.

