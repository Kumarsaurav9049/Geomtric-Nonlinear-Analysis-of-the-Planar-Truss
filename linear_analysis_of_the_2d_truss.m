function linear_analysis_of_the_2d_truss()
example = input(['Enter or Choose an option (Homework question 1 and 2 for standard example solved' ...
            ',\nExample 3.2 of the book David Hutton (Fundamental of finite elemnt analysis)' ...
            ' \nfor custom data press any numeric value other than the (1,2,3)): ']);
        fprintf('solving the standard example: %f \n',example)
        if example == 1
                % Establish coordiantes and global dof's of nodes,
                % node = [name, dof x, x, dof y, y]
                node = [ 1,  1,  1.5,   2,   0 ;
                    2,  3,  3.5,   4,   0 ;
                    3,  5,    0,   6,   6 ;
                    4,  7,    6,   8,   6 ];
                n_node = size(node,1) ;  % number of nodes
                % Discretization - assign elements to global nodes,
                % elem = [name, node i, node j]
                elem = [ 1, node(3,:), node(4,:) ;
                    2, node(1,:), node(4,:) ;
                    3, node(2,:), node(4,:) ];
                n_elem = size(elem,1) ;  % number of elements
                % Geometric and material properties:
                E(1:n_elem) = 200e6;  % Young's modulus (N/m^2)
                A(1:n_elem) = 3480e-6 ;  % Cross-sectional area (m^2)
    
                % Establish applied loads:
                P = sym('P',[2*n_node 1]);   % Load vector as symbolic array
                P(7) =    0;   % load applied to dof 3 (N)
                P(8) = -480;   % load applied to dof 4 (N)
    
                % Apply displacement boundary conditions:
                u = sym('u',[2*n_node 1]);   % Displacement vector as symbolic array
                u(1:4,1) = 0;   % Pinned nodes are 1 and 2
                u(5:6,1) = 0;   % Pinned node 3
        elseif example==2
                node = [ 1, 1,   0, 2,   0 ;
                         2, 3,   3, 4,   0 ;
                         3, 5,   3, 6,   4 ;
                         4, 7,   0, 8,   4 ];
                n_node = size(node,1);   % number of nodes
                % Discretization - assign elements to global nodes,
                % elem = [name, node i, node j]
                elem = [ 1, node(1,:), node(2,:) ;
                         2, node(1,:), node(3,:) ;
                         3, node(4,:), node(2,:) ;
                         4, node(4,:), node(3,:) ;
                         5, node(2,:), node(3,:)];
                n_elem = size(elem,1);  % number of elements
    
                % Geometric and material properties:
                E(1:5) = 200e6;  % Young's modulus (N/m^2)
                A(1:5) = 100e-4 ;  % Cross-sectional area (m^2)
    
                % Establish applied loads:
                P = sym('P',[2*n_node 1]);   % Load vector as symbolic array
                P(3) =    0;   % load applied to dof 3 (N)
                P(4) =    0;   % load applied to dof 4 (N)
                P(5) =    0;   % load applied to dof 5 (N)
                P(6) =  -100;   % load applied to dof 6 (N)
    
                % Apply displacement boundary conditions:
                u = sym('u',[2*n_node 1]);   % Displacement vector as symbolic array
                u(1:2,1) = 0;   % dofs 1 to 2 pinned
                u(7:8,1) = 0;   % dofs 7 to 8 pinned
        elseif example==3
            node = [ 1 1 0  2  0;
                2 3 0  4 40;
                3 5 40 6 40];
            n_node = size(node,1);
            elem =  [1 node(1,:) node(3,:);
                2 node(2,:) node(3,:)];
            n_elem = size(elem,1);

            E = [10*10^6;10*10^6];
            A = [1.5;1.5];
            % Establish applied loads:
            P = sym('P',[2*n_node 1]);   % Load vector as symbolic array
            P(5) =    500;   % load applied to dof 5 (N)
            P(6) =    300;   % load applied to dof 6 (N)
            % Apply displacement boundary conditions:
            u = sym('u',[2*n_node 1])   % Displacement vector as symbolic array
            u(1:4,1) = 0   % node 1 and 2 are pinned supports
           elseif example==4
             node = [1 1 0 2 0; 
                               2 3 1 4 2;
                               3 5 2 6 0;
                               4 7 3 8 2;
                               5 9 4 10 0;
                               6 11 5 12 2;
                               7 13 6 14 0;
                               8 15 7 16 2;
                               9 17 8 18 0];
            n_node = size(node,1);   % number of nodes

            % Discretization - assign elements to global nodes,
            % elem = [name, node i, node j]
            disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)] \n');
            % elem = input('First Enter the name of the member thaen the enter the begigning node and the end node of the memeber :\n');
            elem = [1, node(1,:),node(2,:) ;
                        2, node(2,:),node(3,:) ; 
                        3, node(1,:),node(3,:) ;
                        4, node(2,:),node(4,:) ; 
                        5, node(3,:),node(4,:) ; 
                        6, node(4,:),node(5,:) ; 
                        7, node(5,:),node(3,:) ;
                        8, node(4,:),node(6,:) ; 
                        9, node(5,:),node(6,:) ; 
                        10, node(6,:),node(7,:) ; 
                        11, node(7,:),node(5,:) ; 
                        12, node(6,:),node(8,:) ; 
                        13, node(7,:),node(8,:) ; 
                        14, node(8,:),node(9,:) ; 
                        15, node(7,:),node(9,:)];
            n_elem = size(elem,1) ;  % number of elements

            % Geometric and material properties:
           E(1:n_elem) = 30*10^6;  % Young's modulus (N/m^2)
                A(1:n_elem) = [0.02 ;0.02;0.045;0.045;0.02;0.02;0.045;0.045;0.02;0.02;0.045;0.045;0.02;0.02;0.045] ;  % Cross-sectional area (m^2)

            % Establish applied loads:
            P =   [0;0;15;0;0;-5;0;-7;0;0;0;0;0;-10;0;0;0;0];
             % Apply displacement boundary conditions:
             u = sym('u',[2*n_node 1]);   % Displacement vector as symbolic array
             u = sym('u',[2*n_node 1]);   % Displacement vector as symbolic array
            u([1:2, 18]) = 0;   % Setting specific indices to zero
            % Apply displacement boundary conditions:
            % [u, support_conditions] = specify_support_conditions(n_node);
        else
            % Establish coordiantes and global dof's of nodes,
            % node = [name, dof x, x, dof y, y]
            disp('Format of the node data is like [Node_Number, Dof_x, x_coor, Dof_y, y_coor] \n')
            node = input('Enter the value of the node data with all assigned degree of freedom: \n');
            n_node = size(node,1);   % number of nodes

            % Discretization - assign elements to global nodes,
            % elem = [name, node i, node j]
            disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)] \n');
            elem = input('First Enter the name of the member thaen the enter the begigning node and the end node of the memeber :\n');
            n_elem = size(elem,1) ;  % number of elements

            % Geometric and material properties:
            E = input('Enter the Youngs Modulus of Elasticity of each member in the form of the column vector:\n ');  % Young's modulus (KN/m^2)
            A = input('Enter the Cross-sectional Area of each member in the form of the column vector:\n ');  % Cross-sectional area (m^2)

            % Establish applied loads:
            P = input('Enter the applied Force at the each node in the form of the column vector:\n ');

            % Apply displacement boundary conditions:
            [u, support_conditions] = specify_support_conditions(n_node);
        end
        disp(['*****************************************************' ...
            '********* Performing linear truss analysis...*************' ...
            '****************************']);
        disp(['*****************************************************' ...
            '********* Numerical Solutions...*************' ...
            '*****************************************']);
        % Calculate length and the factors cos and sin for each element:
        for i = 1:n_elem
            L(i) = sqrt( ( elem(i,9) - elem(i,4) )^2 + ( elem(i,11) - elem(i,6) )^2);   % Length of each Member
            c(i) = ( elem(i,9) - elem(i,4) ) / L(i)  ; % Cos(theta)
            s(i) = ( elem(i,11) - elem(i,6) ) / L(i)  ;% Sin(theta)
        end
        % Construct global stiffness matrix from element local matrices:
        K = zeros(2*n_node,2*n_node); % Global stiffness matrix as 2D array of zero's
        for i = 1:n_elem
            k = E(i) * A(i) / L(i) * [     c(i)^2  c(i)*s(i)    -c(i)^2 -c(i)*s(i) ;
                c(i)*s(i)     s(i)^2 -c(i)*s(i)    -s(i)^2 ;
                -c(i)^2 -c(i)*s(i)     c(i)^2  c(i)*s(i) ;
                -c(i)*s(i)    -s(i)^2  c(i)*s(i)     s(i)^2 ]   ;% element local stiffness matrix (kN/mm)

            index = [ elem(i,3) elem(i,5) elem(i,8) elem(i,10) ]   ;% global indices (dofs) for local stiffness matrix
            K(index,index) = K(index,index) + k  ;% add local stiffness matrices to respective global indices (kN/mm)
        end
        % Chop off and find which i-th rows and i-th columns will be kept for the subsequent calculation of 'u':
        not_eliminated = false(2*n_node);    % initiate a false (0) logical array corresponding to the size of K
        for i = 1:2*n_node
            if u(i) ~= 0  % condition for which a row or column is NOT eliminated, i.e., kept
                not_eliminated(i) = true ;   % set the i-th logical entry to true (1)
            end
        end
        u;
        K;
        % Solve for unknown displacements by using the elimination approach with 'not_eliminated' rows and columns:
        u(not_eliminated) = inv(K(not_eliminated,not_eliminated)) * P(not_eliminated) ; % displacement solution (m)
        u = double(u);  % convert symbolic array to double precision
        % Solve for unknown reaction components of the vector 'P':
        P = K * u;  % includes applied loads and reactions (N)
        P = double(P);  % convert symbolic array to double precision
        plot_deformation_of_truss(node, elem, n_node, n_elem, u)
        % plot_deformation_of_truss(node, elem, n_node, n_elem, u)
        [mem_force, strain, stress] = mfss_c(elem, u, E, A, n_elem);
        plot_deformation_and_stresses(node, elem, n_node, n_elem, u, E, A, L);
end