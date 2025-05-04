function plot_deformation_and_stresses(node, elem, n_node, n_elem, u, E, A, L)
    % Call the fss_calc function to get stress, strain, and member force
    [mem_force, strain, stress] = fss_calc(elem, u, E, A, n_elem);

    % Plotting nodes and elements (undeformed and deformed) for stress
    figure(3);
    hold on;
    title('Deformed Truss Structure with Stress Visualization', 'Interpreter', 'latex');

    for i = 1:size(node, 1)
        plot(node(i, 3), node(i, 5), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        text(node(i, 3), node(i, 5), num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
    end

    for i = 1:n_elem
        x1 = node(elem(i, 2), 3);
        y1 = node(elem(i, 2), 5);
        x2 = node(elem(i, 7), 3);
        y2 = node(elem(i, 7), 5);
        mid_x = (x1 + x2) / 2;
        mid_y = (y1 + y2) / 2;
        
        % Determine min and max stresses for coloring elements
        min_sigma = min(stress);
        max_sigma = max(stress);

        % Define a colormap for stress visualization
        colormap('jet'); % You can choose any other colormap as well
        c = colormap;
        
        element_color = c(round(interp1([min_sigma, max_sigma], [1, size(c, 1)], stress(i))), :);
        x1 = node(elem(i, 2), 3) + u(2 * elem(i, 2) - 1);
        y1 = node(elem(i, 2), 5) + u(2 * elem(i, 2));
        x2 = node(elem(i, 7), 3) + u(2 * elem(i, 7) - 1);
        y2 = node(elem(i, 7), 5) + u(2 * elem(i, 7));
        plot([x1, x2], [y1, y2], 'LineWidth', 2, 'Color', element_color);

        % Store legend entry for stress with magnitude
        stress_legend{i} = ['Member ' num2str(elem(i, 1)) ' Stress: ' num2str(stress(i))];
    end

    % Create a legend for stress with all member entries
    legend(stress_legend, 'Location', 'Best', 'Interpreter', 'latex');
    colorbar; % Add a colorbar to show stress scale
    axis on;
    grid on;
    xlabel('X-axis');
    ylabel('Y-axis');

    % Plotting nodes and elements (undeformed and deformed) for member
    % force
    figure(4);
    hold on;
    title('Deformed Truss Structure with Strain Visualization', 'Interpreter', 'latex');

    for i = 1:size(node, 1)
        plot(node(i, 3), node(i, 5), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        text(node(i, 3), node(i, 5), num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
    end

    for i = 1:n_elem
        x1 = node(elem(i, 2), 3);
        y1 = node(elem(i, 2), 5);
        x2 = node(elem(i, 7), 3);
        y2 = node(elem(i, 7), 5);
        mid_x = (x1 + x2) / 2;
        mid_y = (y1 + y2) / 2;
        
        % Determine min and max strains for coloring elements
        min_epsilon = min(strain);
        max_epsilon = max(strain);

        % Define a colormap for strain visualization
        colormap('jet'); % You can choose any other colormap as well
        c = colormap;
        
        element_color = c(round(interp1([min_epsilon, max_epsilon], [1, size(c, 1)], strain(i))), :);
        x1 = node(elem(i, 2), 3) + u(2 * elem(i, 2) - 1);
        y1 = node(elem(i, 2), 5) + u(2 * elem(i, 2));
        x2 = node(elem(i, 7), 3) + u(2 * elem(i, 7) - 1);
        y2 = node(elem(i, 7), 5) + u(2 * elem(i, 7));
        plot([x1, x2], [y1, y2], 'LineWidth', 2, 'Color', element_color);

        % Store legend entry for strain with magnitude
        strain_legend{i} = ['Member ' num2str(elem(i, 1)) ' Strain: ' num2str(strain(i))];
    end

    % Create a legend for strain with all member entries
    legend(strain_legend, 'Location', 'Best', 'Interpreter', 'latex');
    colorbar; % Add a colorbar to show strain scale
    axis on;
    grid on;
    xlabel('X-axis');
    ylabel('Y-axis');

   % Member force plot
figure(5);
hold on;
title('Member force in each element of the truss', 'Interpreter', 'latex');

% Plotting nodes and elements (undeformed and deformed) for member force
for i = 1:size(node, 1)
    plot(node(i, 3), node(i, 5), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(node(i, 3), node(i, 5), num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
end

for i = 1:n_elem
    x1 = node(elem(i, 2), 3);
    y1 = node(elem(i, 2), 5);
    x2 = node(elem(i, 7), 3);
    y2 = node(elem(i, 7), 5);
    mid_x = (x1 + x2) / 2;
    mid_y = (y1 + y2) / 2;
    
    % Determine min and max member forces for coloring elements
    min_mem_force = min(mem_force);
    max_mem_force = max(mem_force);

    % Define a colormap for member force visualization
    colormap('jet'); % You can choose any other colormap as well
    c = colormap;
    
    element_color = c(round(interp1([min_mem_force, max_mem_force], [1, size(c, 1)], mem_force(i))), :);
    x1 = node(elem(i, 2), 3) + u(2 * elem(i, 2) - 1);
    y1 = node(elem(i, 2), 5) + u(2 * elem(i, 2));
    x2 = node(elem(i, 7), 3) + u(2 * elem(i, 7) - 1);
    y2 = node(elem(i, 7), 5) + u(2 * elem(i, 7));
    plot([x1, x2], [y1, y2], 'LineWidth', 2, 'Color', element_color);

    % Display member force magnitude in the middle of the member
    text(mid_x, mid_y, num2str(mem_force(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

    % Store legend entry for member force with magnitude
    mem_force_legend{i} = ['Member ' num2str(elem(i, 1)) ' Force: ' num2str(mem_force(i))];
end

% Create a legend for member force with all member entries
legend(mem_force_legend, 'Location', 'Best', 'Interpreter', 'latex');
colorbar; % Add a colorbar to show member force scale
axis on;
grid on;
xlabel('X-axis');
ylabel('Y-axis');

end