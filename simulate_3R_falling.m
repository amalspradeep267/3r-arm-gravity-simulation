%Simulation of 3R arm falling under gravity under no friction or torque

function simulate_3R_falling()
    %Setting up parameters
    g = 9.81;
    dt = 0.002;         %timestep
    t_final = 10;        %total simulation time
    t = 0:dt:t_final;
    N = length(t);
    
    %Setting Initial conditions
    q = [0; 0; 0]; %theta=q
    q_d = [0; 0; 0]; %angular velocity
    
    %Allocating storage for theta: Is useful in computation
    q_hist = zeros(3, N);
    q_d_hist = zeros(3, N);
    q_dd_hist = zeros(3, N);

    
    %Simulation using Euler integration
    for i = 1:N
        %Store current state
        q_hist(:,i) = q;
        q_d_hist(:,i) = q_d;
        
        %Defining
        M = mass_matrix(q);  %Mass Matrix
        grav = gravity_vector(q, g);  %gravity term
        cor = coriolis_vector(q, q_d);  %Coriolis term
        
        %Forward dynamics: [ M * q_dd = -c - g], for torque=0
        q_dd = M\(-cor - grav); %Solves for theta in terms of above equation
        q_dd_hist(:,i) = q_dd; %Saves value of theta obtained to history
        
        %Euler integration
        q = q + q_d * dt; %gives us the next value for position
        q_d = q_d + q_dd * dt;  %gives us next value for angular velocity
    end

    % Plot results
    plot_results(t, q_hist, q_d_hist, q_dd_hist);
    
    % Animation
    animate_arm(t, q_hist, dt);
end

%Mass matrix
function M = mass_matrix(q)
    q1 = q(1); q2 = q(2); q3 = q(3);

    % Link lengths and masses where all equals 1
    L1=1; L2=1; L3=1;
    m1=1; m2=1; m3=1;

    %Jacobian for mass at end of link 1
    J1 = [ -L1*sin(q1), 0, 0;
            L1*cos(q1), 0, 0];

    %Jacobian for mass at end of link 2
    J2 = [ -L1*sin(q1)-L2*sin(q1+q2), -L2*sin(q1+q2), 0;
            L1*cos(q1)+L2*cos(q1+q2), L2*cos(q1+q2), 0];

    %Jacobian for mass at end of link 3
    J3 = [ -L1*sin(q1)-L2*sin(q1+q2)-L3*sin(q1+q2+q3), -L2*sin(q1+q2)-L3*sin(q1+q2+q3), -L3*sin(q1+q2+q3);
            L1*cos(q1)+L2*cos(q1+q2)+L3*cos(q1+q2+q3), L2*cos(q1+q2)+L3*cos(q1+q2+q3), L3*cos(q1+q2+q3) ];

    %Mass matrix final
    M = m1*(J1'*J1) + m2*(J2'*J2) + m3*(J3'*J3);

    M = M + 1e-4 * eye(3);
end


%% Gravity Vector g(q)
function grav = gravity_vector(q, g)

    q1 = q(1); q2 = q(2); q3 = q(3);
    L1=1; L2=1; L3=1;
    % Calculate the gravity vector, mi=Li=1
    g1 = g * (3* L1 * cos(q1) + 2* L2 * cos(q1 + q2) + L3 * cos(q1 + q2 + q3));
    g2 = g * (2* L2 * cos(q1 + q2) + L3 * cos(q1 + q2 + q3));
    g3 = g * L3 * cos(q1 + q2 + q3);

    %g1 = g * (3*c1 + 2*c12 + c123);
    %g2 = g * (2*c12 + c123);
    %g3 = g * c123;
    
    grav = [g1; g2; g3];
end

%% Coriolis/Centripetal Vector
function cor = coriolis_vector(q, q_d)

    eps = 1e-8;  %For numerical diffrentiation: perturbation
    n = 3;

    %This shows how mass matrix changes with joint angle.
    %i.e partial differentiation of M with joint angles
    dMdq = zeros(n, n, n);  
    
    for i = 1:n
        q_plus = q;  
        q_minus = q;
        q_plus(i)  = q(i) + eps;
        q_minus(i) = q(i) - eps;

        M_plus  = mass_matrix(q_plus);
        M_minus = mass_matrix(q_minus);

        dMdq(:,:,i) = (M_plus - M_minus) / (2*eps);
    end

    %Christoffel-based Coriolis vector
    cor = zeros(n, 1);
    
    for k = 1:n
        for i = 1:n
            for j = 1:n
                Gamma_ijk = 1/2 * ( dMdq(k,j,i) + dMdq(k,i,j) - dMdq(i,j,k));
                cor(k) = cor(k) + Gamma_ijk * q_d(i)* q_d(j);
            end
        end
    end
end

%% Forward Kinematics
function [x, y] = forward_kinematics(q)
    L1 = 1; L2 = 1; L3 = 1;
    
    % Base/joint 1
    x0 = 0;
    y0 = 0;

    
    % Joint 2
    x1 = L1*cos(q(1));
    y1 = L1*sin(q(1));
    
    % Joint 3
    x2 = x1 + L2*cos(q(1) + q(2));
    y2 = y1 + L2*sin(q(1) + q(2));
    
    % End-effector
    x3 = x2 + L3*cos(q(1) + q(2) + q(3));
    y3 = y2 + L3*sin(q(1) + q(2) + q(3));
    
    x = [x0, x1, x2, x3];
    y = [y0, y1, y2, y3];
end

%% Plot Results
function plot_results(t, q_hist, q_d_hist, q_dd_hist)
    figure('Name', '3R Robot Simulation Results', 'Position', [100 100 1200 800]);
    
    % Joint angles
    
    plot(t, q_hist(1,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(t, q_hist(2,:), 'r-', 'LineWidth', 1.5);
    plot(t, q_hist(3,:), 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Angle (rad)');
    title('Joint Angles vs Time');
    legend('\theta_1', '\theta_2', '\theta_3', 'Location', 'best');
    grid on;
    
    figure('Name', '3R Robot Simulation Results', 'Position', [100 100 1200 800]);
    % Joint velocities
    
    plot(t, q_d_hist(1,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(t, q_d_hist(2,:), 'r-', 'LineWidth', 1.5);
    plot(t, q_d_hist(3,:), 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    title('Joint Velocities vs Time');
    legend('$\dot{\theta}_1$', '$\dot{\theta}_2$', '$\dot{\theta}_3$', ...
           'Interpreter', 'latex', 'Location', 'best');
    grid on;
    
    figure('Name', '3R Robot Simulation Results', 'Position', [100 100 1200 800]);
    % Joint accelerations
    
    plot(t, q_dd_hist(1,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(t, q_dd_hist(2,:), 'r-', 'LineWidth', 1.5);
    plot(t, q_dd_hist(3,:), 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Angular Acceleration (rad/s^2)');
    title('Joint Accelerations vs Time');
    legend('$\ddot{\theta}_1$', '$\ddot{\theta}_2$', '$\ddot{\theta}_3$', ...
           'Interpreter', 'latex', 'Location', 'best');
    grid on;
end

%% Animation
function animate_arm(t, q_hist, dt)
    figure('Name', '3R Robot Animation', 'Position', [100 100 800 600]);
    
    % Animation settings for faster simulation
    skip = 50;
    
    for i = 1:skip:length(t)
        [x, y] = forward_kinematics(q_hist(:,i));
        
        clf;
        hold on;
        
        % Drawing links
        plot(x, y, 'b-', 'LineWidth', 3);
        
        % Drawing joints
        plot(x(1), y(1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k');  % Base
        plot(x(2), y(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Joint 2
        plot(x(3), y(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Joint 3
        plot(x(4), y(4), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');  % End-effector
        
        % Drawing masses located at link tip
        plot(x(2), y(2), 'bs', 'MarkerSize', 15, 'MarkerFaceColor', 'b');
        plot(x(3), y(3), 'bs', 'MarkerSize', 15, 'MarkerFaceColor', 'b');
        plot(x(4), y(4), 'bs', 'MarkerSize', 15, 'MarkerFaceColor', 'b');
        
        axis equal;
        xlim([-4 4]);
        ylim([-4 4]);
        grid on;
        xlabel('x (m)');
        ylabel('y (m)');
        title(sprintf('3R Robot Falling Under Gravity (t = %.2f s)', t(i)));
        quiver(3, 3, 0, -0.5, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        text(3.2, 2.7, 'g', 'FontSize', 14);
        
        drawnow;
        pause(dt * skip * 0.2);
    end
end
