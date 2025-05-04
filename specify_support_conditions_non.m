
function [u_non, support_conditions_non] = specify_support_conditions_non(No_non_node)

u_non = sym('u', [2 * No_non_node, 1]);
support_conditions_non = zeros(No_non_node, 1);

disp('Enter support conditions:');
disp('1. Pinned');
disp('2. Roller');
disp('3. None');

for i = 1:No_non_node
    support_choice = input(['Enter support type for Node ' num2str(i) ' (1-Pinned, 2-Roller, 3-None): ']);

    if support_choice == 1
        % Pinned support
        u_non(2 * i - 1:2 * i) = 0;
        support_conditions_non(i) = 1; % Mark node as pinned
    elseif support_choice == 2
        % Roller support
        % For a roller, only the vertical DOF is fixed
        u_non(2 * i) = 0;
        support_conditions_non(i) = 2; % Mark node as roller
    elseif support_choice == 3
        % No support, do nothing
    else
        error('Invalid support type. Please enter 1 (Pinned), 2 (Roller), or 3 (None).');
    end
end

% Display the support conditions for each node
disp('Support Conditions:');
for i = 1:No_non_node
    if support_conditions_non(i) == 1
        disp(['Node ' num2str(i) ' - Pinned']);
    elseif support_conditions_non(i) == 2
        disp(['Node ' num2str(i) ' - Roller']);
    else
        disp(['Node ' num2str(i) ' - None']);
    end
end
end
