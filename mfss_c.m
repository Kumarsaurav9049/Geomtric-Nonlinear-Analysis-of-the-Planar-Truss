function [mem_force, strain, stress] = mfss_c(elem, u, E, A, n_elem)
stress = zeros(n_elem, 1);
strain = zeros(n_elem, 1);
mem_force = zeros(n_elem, 1);


for i = 1:n_elem
    L = sqrt((elem(i, 9) - elem(i, 4))^2 + (elem(i, 11) - elem(i, 6))^2); % Length of each Member
    c = (elem(i, 9) - elem(i, 4)) / L; % Cos(theta)
    s = (elem(i, 11) - elem(i, 6)) / L; % Sin(theta)

    u_1 = u(elem(i, 3), 1) * c + u(elem(i, 5), 1) * s;
    u_2 = u(elem(i, 8), 1) * c + u(elem(i, 10), 1) * s;

    strain(i) = (u_2 - u_1) / L;

    stress(i) = E(i) * strain(i);

    mem_force(i) = stress(i) * A(i);

    % Display results in a long format
    fprintf('Member No: %d\n', i);
    fprintf('Member Length: %.2f\n', L);
    fprintf('Cos(theta): %.4f\n', c);
    fprintf('Sin(theta): %.4f\n', s);
    fprintf('Stress: %.4f\n', stress(i));
    fprintf('Strain: %.4f\n', strain(i));
    fprintf('Member Force: %.4f\n', mem_force(i));
    fprintf('................................\n'); % Add a line break between members

end

end