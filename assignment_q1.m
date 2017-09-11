function assignment_q1
    % Given
    M_o = 15000;      % kg
    M_p = 12000;      % kg
    t_b = 100;           % sec
    u_eq = 3048;      % m/s
    theta_o = 1;        % degree
    delta_t = 0.1;      % sec
    
    % Caculation
    M_b = M_o - M_p; % kg
    massRate = M_p /  t_b; % kg/s
    Thrust = massRate * u_eq;  %N
    R = M_o / M_b;
    
    % Question 1 Solution
    solution('a')
    
    % Question 2 Solution
    solution('b')
    
    % Question 3 Solution
    solution('c')
    
    % Question 4 Solution
    solution('d')
    
    function solution(part)
        D = 0; % N
        g_o = 9.81; % m/s^2
        u_x(1) = 0;  % m/s
        u_y(1) = 0;  % m/s
        u(1) = 0;     % m/s
        u_n(1) = 0; % m/s
        u_r(1) = 0;  % m/s
        x(1) = 0;     % m
        y(1) = 0;     % m
        h(1) = 0;     % m
        theta(1) = theta_o;  % in degree
        M(1) = M_o; % kg
        t(1) = 0;       % s
        for i = 2:1001
            t(i) = t(i - 1) + delta_t;
            if (part == 'a')
                % g = 9.81 (constant)
                g = g_o;
                % Drag = 0
                D = 0;
            elseif (part == 'b')
                % g = 9.81 * (Re/ (Re + h)) 
                g = get_g(g_o, h(i - 1));
                % Drag = 0
                D = 0;
            elseif (part == 'c')
                % g = 9.81 (constant)
                g = g_o;
                % Drag = C_d * (1/2) * rho_h * (u^2) * A_f;
                D = get_D(h(i - 1), u(i -1));
            elseif (part == 'd')
                % g = 9.81 * (Re/ (Re + h)) 
                g = get_g(g_o, h(i - 1));
                 % Drag = C_d * (1/2) * rho_h * (u^2) * A_f;
                D = get_D(h(i - 1), u(i - 1));
            end
            delta_u_tOld = get_delta_u_tOld(massRate, u_eq, M(i -1), D, g, theta(i - 1), delta_t);
            delta_un_tOld = get_delta_un_tOld(g, theta(i - 1), delta_t);
            delta_ur_tOld = get_delta_ur_tOld(delta_u_tOld, delta_un_tOld);
            delta_theta_tOld = get_delta_theta_tOld(delta_un_tOld, u_r(i - 1), delta_u_tOld);
            delta_phi_tOld = get_delta_phi_tOld(delta_un_tOld, delta_u_tOld);
            delta_ux_tOld = get_delta_ux_tOld(delta_ur_tOld, theta(i -1), delta_phi_tOld);
            delta_uy_tOld = get_delta_uy_tOld(delta_ur_tOld, theta(i -1), delta_phi_tOld);

            theta(i) = get_theta_tNew(theta(i - 1), delta_theta_tOld);
            u_x(i) = get_ux_tNew(u_x(i - 1), delta_ux_tOld);
            u_y(i) = get_uy_tNew(u_y(i - 1), delta_uy_tOld);
            u_r(i) = get_ur_tNew(u_x(i), u_y(i));
            u(i) = get_u_tNew(u_r(i));

            delta_x_tOld = get_delta_x_tOld(u_x(i), delta_t);
            delta_y_tOld = get_delta_y_tOld(u_y(i), delta_t);
            delta_M_tOld = get_delta_M_tOld(massRate, delta_t);

            x(i) = get_x_tNew(x(i - 1), delta_x_tOld);
            y(i) = get_y_tNew(y(i - 1), delta_y_tOld);
            M(i) = get_M_tNew(M(i - 1), delta_M_tOld);
            h(i) = get_h_tNew(y(i));
        end
        
        %%%%%%%OUTPUT%%%%%%%%
        
        %%%Answers%%%
        % Burnout height for respective question
        disp(['For Question 1 part ', part])
        disp(['Burnout height(h_b) is ', num2str(h(end)), ' m']);
        % Burnout speed for respective question
        disp(['Burnout speed(v_b) is ', num2str(u(end)), ' m/s']);
        % Burnout theta for respective question
        disp(['Burnout angle(theta) is ', num2str(theta(end)), ' degree']);
        
        %%%Plots%%%
        %%%Plot y-coordinate vs x-coordinate%%%
        save_fig = figure
        hold on
        plot(x, y)
        title(['x vs y for part ', part])
        ylabel('x co-ordinate (in m)')
        xlabel('y co-ordinate (in m)')
        hold off
        saveas(save_fig, ['x_vs_y_part_', part, '.png'])
        %%%Plot speed vs time%%%
        save_fig = figure
        hold on
        plot(t, u)
        title(['u vs t for part ', part])
        ylabel('Velocity (in m/s)')
        xlabel('Time (in sec)')
        hold off
        saveas(save_fig, ['u_vs_t_part_', part, '.png'])
        %%%Plot elevation angle vs time%%%
        save_fig = figure
        hold on
        plot(t, theta)
        title(['theta vs t for part ', part])
        ylabel('Elevation angle (in degree)')
        xlabel('Time (in sec)')
        hold off
        saveas(save_fig, ['theta_vs_t_part_', part, '.png'])
        %%%Plot height vs time%%
        save_fig = figure
        hold on
        plot(t, h)
        title(['height vs t for part ', part])
        ylabel('Height (in m)')
        xlabel('Time (in sec)')
        hold off
        saveas(save_fig, ['height_vs_t_part_', part, '.png'])
    end
     
    %%%Functions%%%
    function g = get_g(g_o, h)
        g = g_o * (6400000 / (6400000 + h));
    end

    function D = get_D(h, u)
        rho_h = 1.2 * exp(-2.9 * (10^(-5)) * (h^(1.15)));
        A_f = 1;
        C_d = 0.1;
        D = C_d * (1/2) * rho_h * (u^2) * A_f;
    end

    function delta_u_tOld = get_delta_u_tOld(massRate_tOld, ueq_tOld, M_tOld, D_tOld, g_tOld, theta_tOld, delta_t)
        delta_u_tOld = ( ((massRate_tOld * ueq_tOld) / M_tOld) - (D_tOld / M_tOld) - (g_tOld * cosd(theta_tOld)) ) * delta_t;
    end

    function delta_un_tOld = get_delta_un_tOld(g_tOld, theta_tOld, delta_t)
        delta_un_tOld = g_tOld * sind(theta_tOld) * delta_t;
    end

    function delta_ur_tOld = get_delta_ur_tOld(delta_u_tOld, delta_un_tOld)
        delta_ur_tOld = ( (delta_u_tOld)^2 + (delta_un_tOld)^2) ^ (1/2);
    end

    function delta_theta_tOld = get_delta_theta_tOld(delta_un_tOld, ur_tOld, delta_u_tOld)
        delta_theta_tOld = atand((delta_un_tOld) / (ur_tOld + delta_u_tOld));
    end

    function delta_phi_tOld = get_delta_phi_tOld(delta_un_tOld, delta_u_tOld)
        delta_phi_tOld = atand(delta_un_tOld / delta_u_tOld);
    end
    
    function delta_ux_tOld = get_delta_ux_tOld(delta_ur_tOld, theta_tOld, delta_phi_tOld)
        delta_ux_tOld =delta_ur_tOld * (sind(theta_tOld + delta_phi_tOld));
    end
    
   function delta_uy_tOld = get_delta_uy_tOld(delta_ur_tOld, theta_tOld, delta_phi_tOld)
        delta_uy_tOld =delta_ur_tOld * (cosd(theta_tOld + delta_phi_tOld));
   end

    function theta_tNew = get_theta_tNew(theta_tOld, delta_theta_tOld)
        theta_tNew = theta_tOld + delta_theta_tOld;
    end

    function  ux_tNew = get_ux_tNew(ux_tOld, delta_ux_tOld)
        ux_tNew = ux_tOld + delta_ux_tOld;
    end

    function  uy_tNew = get_uy_tNew(uy_tOld, delta_uy_tOld)
        uy_tNew = uy_tOld + delta_uy_tOld;
    end

    function ur_tNew = get_ur_tNew(ux_tNew, uy_tNew)
        ur_tNew = ((ux_tNew ^ 2) + (uy_tNew ^ 2))^(1/2);
    end

    function u_tNew = get_u_tNew(ur_tNew)
        u_tNew = ur_tNew;
    end

    function delta_x_tOld = get_delta_x_tOld(ux_tNew, delta_t)
        delta_x_tOld = ux_tNew * delta_t;
    end

    function delta_y_tOld = get_delta_y_tOld(uy_tNew, delta_t)
        delta_y_tOld = uy_tNew * delta_t;
    end

    function delta_M_tOld = get_delta_M_tOld(massRate_tOld, delta_t)
        delta_M_tOld = -massRate_tOld * delta_t;
    end

    function x_tNew = get_x_tNew(x_tOld, delta_x_tOld)
        x_tNew = x_tOld + delta_x_tOld;
    end
    
   function y_tNew = get_y_tNew(y_tOld, delta_y_tOld)
        y_tNew = y_tOld + delta_y_tOld;
   end

    function M_tNew = get_M_tNew(M_tOld, delta_M_tOld)
        M_tNew = M_tOld + delta_M_tOld;
    end

    function h_tNew = get_h_tNew(y_tNew)
        h_tNew = y_tNew;
    end

end
