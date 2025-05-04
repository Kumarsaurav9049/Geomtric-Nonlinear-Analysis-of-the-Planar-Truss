clc
clear all
close all
node = [ 1,  4,    0,   5,   0 ;
         2,  1,    4,   2,   3 ;
         3,  3,    8,   6,   0  ];    %input('Enter the value of the node data with all assigned degree of freedom: \n')
n_node = size(node,1);   % number of nodes

elem = [ 1, node(1,:), node(2,:) ;
         2, node(3,:), node(2,:) ;
         3, node(1,:), node(3,:) ];
n_elem = size(elem,1) ;  % number of elements
E = [70*10^6;70*10^6;70*10^6]; %input('Enter the Youngs Modulus of Elasticity of each member in the form of the column vector:\n ');  % Young's modulus (KN/m^2)
A = [645*10^-6;645*10^-6;645*10^-6]; %input('Enter the Cross-sectional Area of each member in the form of the column vector:\n ');  % Cross-sectional area (m^2)

       

for i = 1:n_elem
    L(i) = sqrt( ( elem(i,9) - elem(i,4) )^2 + ( elem(i,11) - elem(i,6) )^2) ;  % Length of each Member
    c(i) = ( elem(i,9) - elem(i,4) ) / L(i);   % Cos(theta)
    s(i) = ( elem(i,11) - elem(i,6) ) / L(i);  % Sin(theta)
end
non_node = [ 1,  4,    (0),   5,   (0);
    2,  1,    (4 + 0.1181),   2,   (3 -0.4651);
    3,  3,    (8 + 0.2362),   6,   (0 ) ]
No_non_node = size(non_node,1);
non_elem = [ 1, non_node(1,:), non_node(2,:) ;
    2, non_node(3,:), non_node(2,:) ;
    3, non_node(1,:), non_node(3,:) ]
No_non_elem = size(non_elem,1);

local_tangent_stiffness_matrices = zeros(4, 4, No_non_elem);
% Construct global stiffness matrix from element local matrices:
K_t = zeros(2*No_non_node,2*No_non_node) ; % Global stiffness matrix as 2D array of zero's
F_global = zeros(2*No_non_node,1) ;
for i = 1:No_non_elem
    L_prime(i)  = sqrt( ( non_elem(i,9) - non_elem(i,4) )^2 + ( non_elem(i,11) - non_elem(i,6) )^2)  % Length of each Member
    c_x(i)      = ( non_elem(i,9) - non_elem(i,4) ) / L_prime(i)   % Cos(theta)
    c_y(i)      = ( non_elem(i,11) - non_elem(i,6) ) / L_prime(i) % Sin(theta)
    U(i)        = L(i)-L_prime(i)
    Q(i)        = (A(i)*E(i)/L(i))*U(i)
    T           = [c_x(i) c_y(i) -c_x(i) -c_y(i)];
    F           = T'*Q(i)
    k_t         = E(i) * A(i) / L(i) *T'*T  + (Q(i)/L_prime(i))*[ -c_y(i)^2         c_x(i)*c_y(i)    c_y(i)^2         -c_x(i)*c_y(i);
        c_x(i)*c_y(i)     -c_x(i)^2       -c_x(i)*c_y(i)            c_x(i)^2;
        c_y(i)^2         -c_x(i)*c_y(i)   -c_y(i)^2         c_x(i)*c_y(i);
        -c_x(i)*c_y(i)        c_x(i)^2     c_x(i)*c_y(i)     -c_x(i)^2  ]  % element local stiffness matrix (kN/mm)
    local_tangent_stiffness_matrices(:, :, i) = k_t;
    Findex = [elem(i,3); elem(i,5); elem(i,8); elem(i,10) ];
    F_global(Findex) = F_global(Findex) + F
    index = [ elem(i,3) elem(i,5) elem(i,8) elem(i,10) ] ;  % global indices (dofs) for local stiffness matrix
    K_t(index,index) = K_t(index,index) + k_t  % add local stiffness matrices to respective global indices (kN/mm)
end
% for i = 1:n_elem
%     disp("Local Stiffness Matrix for element " + num2str(i) + ":");
%     disp(local_stiffness_matrices(:, :, i))
% end
[u_non, support_conditions_non] = specify_support_conditions_non(No_non_node)
% Chop off and find which i-th rows and i-th columns will be kept for the subsequent calculation of 'u':
not_eliminated = false(2*3) ;   % initiate a false (0) logical array corresponding to the size of K
for i = 1:2*No_non_node
    if u_non(i) ~= 0  % condition for which a row or column is NOT eliminated, i.e., kept
        not_eliminated(i) = true    % set the i-th logical entry to true (1)
    end
end
u_non;
K_t;
% f1  = [0;0;0;-1554.98;70.322;0] ;
Net_force = [0;0;0;-2000;0;0] - F_global
% Solve for unknown displacements by using the elimination approach with 'not_eliminated' rows and columns:
u_non(not_eliminated) = inv(K_t(not_eliminated,not_eliminated)) * Net_force(not_eliminated); % displacement solution (m)
u_non = double(u_non)% convert symbolic array to double precision
% utemp = u_non
% total_displacement = u + u_non
convergence_criterion = sqrt(sum(u_non.^2)/(sum(u.^2)))
% end