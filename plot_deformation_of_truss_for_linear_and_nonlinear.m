function plot_deformation_of_truss_for_linear_and_nonlinear(node, elem, n_node, n_elem, u_linear,u_nonlinear)
% Plotting nodes (undeformed)
figure;
hold on;
title(['Undeformed and Deformed Truss Structure [Dashed Blue' ...
    ' Line Repersents the undeformed truss and red repersents the deformed truss ' ...
    '(On true scale) ]'],Interpreter='latex');

for i = 1:size(node, 1)
    plot(node(i, 3), node(i, 5),'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(node(i, 3), node(i, 5), num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
    % undeformed_node_legend{i} = ['undeformed nodes' num2str(node(i, 1))];

end
% legend(undeformed_node_legend,'Location','best')
% Plotting elements (undeformed)
for i = 1:size(elem, 1)
    x1 = node(elem(i, 2), 3);
    y1 = node(elem(i, 2), 5);
    x2 = node(elem(i, 7), 3);
    y2 = node(elem(i, 7), 5);
    plot([x1, x2], [y1, y2], '--b', 'LineWidth', 2,'HandleVisibility','off');
    mid_x = (x1 + x2) / 2;
    mid_y = (y1 + y2) / 2;
    text(mid_x, mid_y, num2str(elem(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 20);
end
% Plotting nodes (deformed)
for i = 1:n_node
    deformed_x = node(i, 3) + u_linear(2 * i - 1);
    deformed_y = node(i, 5) + u_linear(2 * i);
    plot(deformed_x, deformed_y, 'o', 'MarkerSize', 10,'HandleVisibility','off');
    text(deformed_x, deformed_y, num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
end

% Plotting elements (deformed)
for i = 1:n_elem
    x1 = node(elem(i, 2), 3) + u_linear(2 * elem(i, 2) - 1);
    y1 = node(elem(i, 2), 5) + u_linear(2 * elem(i, 2));
    x2 = node(elem(i, 7), 3) + u_linear(2 * elem(i, 7) - 1);
    y2 = node(elem(i, 7), 5) + u_linear(2 * elem(i, 7));
    plot([x1, x2], [y1, y2], 'r', 'LineWidth', 2);
end


% Plotting nodes (Nonlinear deformed)
for i = 1:n_node
    deformed_x = node(i, 3) + u_nonlinear(2 * i - 1);
    deformed_y = node(i, 5) + u_nonlinear(2 * i);
    plot(deformed_x, deformed_y, 'o', 'MarkerSize', 10,'HandleVisibility','off');
    text(deformed_x, deformed_y, num2str(node(i, 1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
end

% Plotting elements (deformed)
for i = 1:n_elem
    x1 = node(elem(i, 2), 3) + u_nonlinear(2 * elem(i, 2) - 1);
    y1 = node(elem(i, 2), 5) + u_nonlinear(2 * elem(i, 2));
    x2 = node(elem(i, 7), 3) + u_nonlinear(2 * elem(i, 7) - 1);
    y2 = node(elem(i, 7), 5) + u_nonlinear(2 * elem(i, 7));
    plot([x1, x2], [y1, y2], 'g', 'LineWidth', 2);
end

axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
legend('Linearly Deformed Truss Members','Linearly Deformed Truss Members','Linearly Deformed Truss Members', ...
    'Geometrically Nonlinearly deformed Truss Members','Geometrically Nonlinearly deformed Truss Members', ...
    'Geometrically Nonlinearly deformed Truss Members',Interpreter='latex')

end
