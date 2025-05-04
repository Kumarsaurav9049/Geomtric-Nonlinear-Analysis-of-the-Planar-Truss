function nonlinear_analysis_of_the_2d_truss()
% Establish coordiantes and global dof's of nodes,
            % node = [name, dof x, x, dof y, y]
            example = input(['Enter or Choose an option (Chapter 10 of the Aslam Kassimals Book Matrix Analysis of Structure is titled'... 
                           '\nGeometrically Nonlinear Analysis of Plane Trusses Problems (10.2) and (10.5) in the book have been solved' ...
                           ' \nas the usual example To view the outcomes of the standard problem, press 1 and 2, respectively none of the' ...
                           ' \nmore significant HOMEWORK QUESTION answered to visualise the home work problem solution press 3 and to solve' ...
                           ' \nany general plane truss, press any value except (1,2,3).:\n'] );
            fprintf('solving the standard example: %f \n',example)
                if example == 1
                        % disp('Format of the node data is like [Node_Number, Dof_x, x_coor, Dof_y, y_coor] ')
                        node = [ 1,  1,    0,   2,   0 ;
                            2,  3,    4,   4,   3 ;
                            3,  5,    8,   6,   0 ];
                        % input('Enter the value of the node data with all assigned degree of freedom: \n')
                        n_node = size(node,1);   % number of nodes
            
                        % Discretization - assign elements to global nodes,
                        % elem = [name, node i, node j]
                        % disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)]')
                        elem = [ 1, node(1,:), node(2,:) ;
                            2, node(1,:), node(3,:) ;
                            3, node(2,:), node(3,:) ];
                        % input('First Enter the name of the member thaen the enter the begigning node and the end node of the memeber :\n')
                        n_elem = size(elem,1) ;  % number of elements
            
                        % Geometric and material properties:
                        E = [70*10^6;70*10^6;70*10^6]; %input('Enter the Youngs Modulus of Elasticity of each member in the form of the column vector:\n ');  % Young's modulus (KN/m^2)
                        A = [645*10^-6;645*10^-6;645*10^-6]; %input('Enter the Cross-sectional Area of each member in the form of the column vector:\n ');  % Cross-sectional area (m^2)
            
                        % Establish applied loads:
                        P = [0;0;0;-2000;0;0]; %input('Enter the applied Force at the each node in the form of the column vector:\n ')
            
                        % Apply displacement boundary conditions:
                        % u = sym('u', [2 * n_node, 1]);
                        % u(1:3) = 0;
                        [u, support_conditions] = specify_support_conditions(n_node);
                    elseif example == 2
            
                        % node = [name, dof x, x, dof y, y]
                        % disp('Format of the node data is like [Node_Number, Dof_x, x_coor, Dof_y, y_coor] ')
                        node = [ 1,  1,     0,    2,    0 ;
                                 2,  3,     4,    4,    0 ;
                                 3,  5,    16,   6,   16  ];    %input('Enter the value of the node data with all assigned degree of freedom: \n')
                        n_node = size(node,1);   % number of nodes
            
                        % Discretization - assign elements to global nodes,
                        % elem = [name, node i, node j]
                        % disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)]')
                        elem = [ 1, node(1,:), node(2,:) ;
                                 2, node(1,:), node(3,:) ;
                                 3, node(2,:), node(3,:) ];
                        % input('First Enter the name of the member thaen the enter the begigning node and the end node of the memeber :\n')
                        n_elem = size(elem,1) ;  % number of elements
            
                        % Geometric and material properties:
                        E = [70*10^6;70*10^6;70*10^6]; % input('Enter the Youngs Modulus of Elasticity of each member in the form of the column vector:\n ');  % Young's modulus (KN/m^2)
                        A = [1200*10^-6;1200*10^-6;1200*10^-6]; %input('Enter the Cross-sectional Area of each member in the form of the column vector:\n ');  % Cross-sectional area (m^2)
            
                        % Establish applied loads:
                        P = [0;0;0;0;0;-150]; %input('Enter the applied Force at the each node in the form of the column vector:\n ')
            
                        % Apply displacement boundary conditions:
                        % u = sym('u', [2 * n_node, 1]);
                        % u(1:3) = 0;
                        [u, support_conditions] = specify_support_conditions(n_node);
                elseif example==3
                        % disp('Format of the node data is like [Node_Number, Dof_x, x_coor, Dof_y, y_coor] ')
                        % node = [name, dof x, x, dof y, y]   
                        yy = input('Enter the Roll number of the Student:- \n');
                        coordinate_xy_value_1 = 1.5*(1 + 0.01*yy);
                        coordinate_xy_value_2 = 2*(1 + 0.01*yy);
                        node = [ 1,  1,                          0,     2,   0 ;
                                 2,  3,    2*coordinate_xy_value_2,     4,   0;
                                 3,  5,      coordinate_xy_value_2,     6,   coordinate_xy_value_1 ;
                                 4,  7,      coordinate_xy_value_2,     8,   2*coordinate_xy_value_1 ];
                        % input('Enter the value of the node data with all assigned degree of freedom: \n')
                        n_node = size(node,1);   % number of nodes
            
                        % Discretization - assign elements to global nodes,
                        % elem = [name, node i, node j]
                        % disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)]')
                        elem = [ 1, node(1,:), node(3,:) ;
                                 2, node(2,:), node(3,:) ;
                                 3, node(4,:), node(3,:) ];
                        % input('First Enter the name of the member thaen the enter the begigning node and the end node of the memeber :\n')
                        n_elem = size(elem,1) ;  % number of elements
            
                        % Geometric and material properties:
                        E = [9*10^6;9*10^6;9*10^6]; %input('Enter the Youngs Modulus of Elasticity of each member in the form of the column vector:\n ');  % Young's modulus (KN/m^2)
                        A = [880*10^-6;880*10^-6;440*10^-6]; %input('Enter the Cross-sectional Area of each member in the form of the column vector:\n ');  % Cross-sectional area (m^2)
            
                        % Establish applied loads:
                        P = [0;0;0;0;1600;-3200;0;0]; %input('Enter the applied Force at the each node in the form of the column vector:\n ')
            
                        % Apply displacement boundary conditions:
                        u = sym('u', [2 * n_node, 1]);
                        u([1:4, 7:8]) = 0;

                        % [u, support_conditions] = specify_support_conditions(n_node);
               
                else
                        disp('Format of the node data is like [Node_Number, Dof_x, x_coor, Dof_y, y_coor] ');
                        node = input('Enter the value of the node data with all assigned degree of freedom: \n');
                        n_node = size(node,1);   % number of nodes
            
                        % Discretization - assign elements to global nodes,
                        % elem = [name, node i, node j]
                        disp('Format of the connectivity matrix for the element data is like [name, beigning_node(i,:), end_node(j,:)]');
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
                tolereance = input('Enter the Value up to What decimal precision result is required try to put maximum decimal point for good result (example : 10e-6) ');

                    disp(['*****************************************************' ...
                        '********* First Performing linear truss analysis...*************' ...
                        '****************************']);
                    disp(['*****************************************************' ...
                        '********* Numerical Solutions...*************' ...
                        '*****************************************']);
                    % Calculate length and the factors cos and sin for each element:
                    for i = 1:n_elem
                        L(i) = sqrt( ( elem(i,9) - elem(i,4) )^2 + ( elem(i,11) - elem(i,6) )^2) ;  % Length of each Member
                        c(i) = ( elem(i,9) - elem(i,4) ) / L(i) ;  % Cos(theta)
                        s(i) = ( elem(i,11) - elem(i,6) ) / L(i);  % Sin(theta)
                    end
            
                    local_stiffness_matrices = zeros(4, 4, n_elem);
            
                    % Construct global stiffness matrix from element local matrices:
                    K = zeros(2*n_node,2*n_node) ; % Global stiffness matrix as 2D array of zero's
                        for i = 1:n_elem
                            k = E(i) * A(i) / L(i) * [     c(i)^2  c(i)*s(i)    -c(i)^2 -c(i)*s(i) ;
                                c(i)*s(i)     s(i)^2 -c(i)*s(i)    -s(i)^2 ;
                                -c(i)^2 -c(i)*s(i)     c(i)^2  c(i)*s(i) ;
                                -c(i)*s(i)    -s(i)^2  c(i)*s(i)     s(i)^2 ]   ;% element local stiffness matrix (kN/mm)
                            local_stiffness_matrices(:, :, i) = k;
                            index = [ elem(i,3) elem(i,5) elem(i,8) elem(i,10) ]   ;% global indices (dofs) for local stiffness matrix
                            K(index,index) = K(index,index) + k  ;% add local stiffness matrices to respective global indices (kN/mm)
                        end
                        for i = 1:n_elem
                            % disp("Local Stiffness Matrix for element " + num2str(i) + ":");
                            % disp(local_stiffness_matrices(:, :, i))
                        end
                        % Chop off and find which i-th rows and i-th columns will be kept for the subsequent calculation of 'u':
                        not_eliminated = false(2*n_node) ;   % initiate a false (0) logical array corresponding to the size of K
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
                    utemp = u;
                    u_linear = utemp;
                    disp('linear Analyis Result for the deflection')
                    disp(utemp)
                    Rf_linear = K*utemp;
                    [~, col] = size(node);
            
                    % Create an empty matrix for non_node
                    non_node = [];
            
                    % Loop through each row of the node matrix
                        for i = 1:size(node, 1)
                            non_node_row = [node(i, 1), node(i, 2), (node(i, 3) + utemp(2*i-1)), node(i, 4), (node(i, 5) + utemp(2*i))];
                            non_node = [non_node; non_node_row];
                        end
                   disp(['*****************************************************' ...
                        '********* Now Performing nonlinear truss analysis...*************' ...
                        '****************************']);
                    disp(['*****************************************************' ...
                        '********* Numerical Solutions...*************' ...
                        '*****************************************']);
                % Display the resulting non_node matrix
                % disp("Non-node matrix:")
                % disp(non_node);
        
                % non_node = [ 1,  1,    (0 + utemp(1,1)),   2,   (0 + utemp(2,1));
                %              2,  3,    (4 + utemp(3,1)),   4,   (3 + utemp(4,1));
                %              3,  5,    (8 + utemp(5,1)),   6,   (0 + utemp(6,1)) ];
                No_non_node = size(non_node,1);
                % [u, ~] = specify_support_conditions_non(No_non_node);
        
                convergence_criterion = inf; % Assuming an initial value
                max_iterations = 100000; % Set the maximum number of iterations
                convergence_data =[];
                iteration = 0; % Initialize the iteration counter
                while convergence_criterion >= tolereance && iteration < max_iterations
                    % Coomon code
                    % non_node = [ 1,  1,    (0 + utemp(1,1)),   2,   (0 + utemp(2,1));
                    %              2,  3,    (4 + utemp(3,1)),   4,   (3 + utemp(4,1));
                    %              3,  5,    (8 + utemp(5,1)),   6,   (0 + utemp(6,1)) ];
                    % Generalized Code
                    [~, col] = size(node);
        
                    % Create an empty matrix for non_node
                    non_node = [];
        
                    % Loop through each row of the node matrix
                    for i = 1:size(node, 1)
                        non_node_row = [node(i, 1), node(i, 2), (node(i, 3) + utemp(2*i-1)), node(i, 4), (node(i, 5) + utemp(2*i))];
                        non_node = [non_node; non_node_row];
                    end
                    No_non_node = size(non_node,1);
                    % Common Code
                    % if example == 1
                    %       non_elem = [   1, non_node(1,:), non_node(2,:) ;
                    %                      2, non_node(1,:), non_node(3,:) ;
                    %                      3, non_node(2,:), non_node(3,:) ];
                    % elseif example==2
                    %       non_elem = [   1, non_node(1,:), non_node(2,:) ;
                    %                      2, non_node(1,:), non_node(3,:) ;
                    %                      3, non_node(2,:), non_node(3,:) ];
                    % elseif example==3
                    % non_elem = [ 1, non_node(1,:), non_node(3,:) ;
                    %              2, non_node(2,:), non_node(3,:) ;
                    %              3, non_node(4,:), non_node(3,:) ];
                    % else
                    % disp('Format of the connectivity matrix in [Elem. no. non_node(1,:) non_node(2,:)')
                    non_elem = [elem(:,1) non_node(elem(:,2),:) non_node(elem(:,7),:)];%input('Eenter the connectivity matrix what you have entered in the above but this time write with the non_node matrix data\n')
                    % end
                    %generalized code for the non_elem
                    No_non_elem = size(non_elem,1);
        
                    local_tangent_stiffness_matrices = zeros(4, 4, n_elem);
        
                    % Construct global stiffness matrix from element local matrices:
        
                    K_t = zeros(2*No_non_node,2*No_non_node) ; % Global stiffness matrix as 2D array of zero's
                    F_g = zeros(2*No_non_node,1) ;
                    for i = 1:No_non_elem
                        L_prime(i)  = sqrt( ( non_elem(i,9) - non_elem(i,4) )^2 + ( non_elem(i,11) - non_elem(i,6) )^2) ; % Length of each Member
                        c_x(i)      = ( non_elem(i,9) - non_elem(i,4) ) / L_prime(i) ;  % Cos(theta)
                        c_y(i)      = ( non_elem(i,11) - non_elem(i,6) ) / L_prime(i); % Sin(theta)
                        U(i)        = L(i)-L_prime(i);
                        Q(i)        = (A(i)*E(i)/L(i))*U(i);
                        T           = [c_x(i) c_y(i) -c_x(i) -c_y(i)];
                        F_local     = T'*Q(i);
                        k_t         = E(i) * A(i) / L(i) *T'*T  + (Q(i)/L_prime(i))*[ -c_y(i)^2         c_x(i)*c_y(i)    c_y(i)^2         -c_x(i)*c_y(i);
                            c_x(i)*c_y(i)     -c_x(i)^2       -c_x(i)*c_y(i)            c_x(i)^2;
                            c_y(i)^2         -c_x(i)*c_y(i)   -c_y(i)^2         c_x(i)*c_y(i);
                            -c_x(i)*c_y(i)        c_x(i)^2     c_x(i)*c_y(i)     -c_x(i)^2  ]  ;% element local stiffness matrix (kN/mm)
                        local_tangent_stiffness_matrices(:, :, i) = k_t;
                        index = [ elem(i,3) elem(i,5) elem(i,8) elem(i,10) ]   ;% global indices (dofs) for local stiffness matrix
                        F_g(index') = F_g(index') + F_local;
                        K_t(index,index) = K_t(index,index) + k_t ; % add local stiffness matrices to respective global indices (kN/mm)
                    end

                    % for i = 1:n_elem
                    %     disp("Local Stiffness Matrix for element " + num2str(i) + ":");
                    %     disp(local_stiffness_matrices(:, :, i))
                    % end
                    % % Chop off and find which i-th rows and i-th columns will be kept for the subsequent calculation of 'u':
                    not_eliminated = false(2*3) ;   % initiate a false (0) logical array corresponding to the size of K
                    for i = 1:2*No_non_node
                        if u(i) ~= 0  % condition for which a row or column is NOT eliminated, i.e., kept
                            not_eliminated(i) = true  ;  % set the i-th logical entry to true (1)
                        end
                    end
                    u;
                    K_t;
                    Net_force = P-F_g; % Residual Force
                    % Solve for unknown displacements by using the elimination approach with 'not_eliminated' rows and columns:
                    u(not_eliminated) = inv(K_t(not_eliminated,not_eliminated)) * Net_force(not_eliminated); % displacement solution (m)
                    u = double(u);% convert symbolic array to double precision
                    convergence_criterion = sqrt(sum(u.^2)/(sum(utemp.^2)));
                    utemp = u + utemp;
                    iteration = iteration + 1;
                    u_nonlinear = utemp;
                    convergence_data = [convergence_data;iteration,convergence_criterion];
                    plot(convergence_data(:,1),convergence_data(:,2),'-o', 'LineWidth', 2);
                    xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman','FontWeight','bold');
                    ylabel('Convergence Criteria', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman','FontWeight','bold');
                    title('Convergence Plot', 'Interpreter', 'latex','FontWeight','bold', 'FontSize', 14, 'FontName', 'Times New Roman')
                    latexExpression = 'Convergence Criteria  $=\sqrt{\frac{\sum\limits_{j=1}^{n} (\Delta d_j)^2}{\sum\limits_{j=1}^{n}(d_j)^2}}$';
                    text('Interpreter', 'latex', 'String', latexExpression, 'Units', 'normalized', 'Position', [1, 1], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontName', 'Times New Roman');            grid on
                end
                disp('The Maximum no of iteratioons for your desired tolerance ');
                disp(iteration)
                disp('The Final Deformation after the  geometric nonlinear analsis');
                disp(u_nonlinear)
                % Final Member Force Calculationin the Truss
                disp('The Member Forces are : Tension (-) and Compressive (+)');
                disp(Q)
                plot_deformation_of_truss_for_linear_and_nonlinear(node, elem, n_node, n_elem, u_linear,u_nonlinear)             
end
