function make_plots_pODE(yns,zns,iter_vec,plots_params)
	% function make_plots_pODE(yns,zns,iter_vec,plots_params)
	%
	% Make plots for partitioned ODEs.
	%
	% Inputs:
	%	yns, zns: The matrices consisting of the update steps (yn,zn)
	%	iter_vec: vector of iterations
	%	plots_params: Structured array with the following fields:
	%		h: stepsize
	%		N: number of steps 
	%		s: number of stages
	%		problem_str: describes the partitioned ODE to be solved.
	%		method_str: describes the PRK method
	%		plot_str: sets a plotting option.
	%		iterative_scheme: string describing the iterative method
	%		plot_iterations: Boolean value determining whether to make 
	%			plot(s) of the iterations.
	%
	% Outputs:
	%	None.

	% unpack plots_params
	h = plots_params.h;
	N = plots_params.N;
	s = plots_params.s;
	problem_str = plots_params.problem_str;
	method_str = plots_params.method_str;
	plot_str = plots_params.plot_str;
	iterative_scheme = plots_params.iterative_scheme;
	plot_iterations = plots_params.plot_iterations; 

	if strcmpi(problem_str, 'pendulum')
		E = plots_params.E; % energy
	elseif strcmpi(problem_str, 'pendulum theta')
		l = plots_params.pendulum_struct.l;
		g0 = plots_params.pendulum_struct.g0;
	elseif strcmpi(problem_str, 'double pendulum angles')
		l1 = plots_params.pendulum_struct.l1;
		l2 = plots_params.pendulum_struct.l2;
		g0 = plots_params.pendulum_struct.g0;
	elseif strcmpi(problem_str, 'two body')
		alpha0 = plots_params.two_body_struct.alpha0;
		T = plots_params.two_body_struct.T;
		a = plots_params.two_body_struct.a;
		e = plots_params.two_body_struct.e;
	end

	ny = size(yns,1);
	nz = size(zns,1);

	% --------------------
	% Plots of iterations 
	% --------------------

	if plot_iterations

		line2 = sprintf(', $s=%i$, $h=%0.3f$, $N=%i$', s, h, N);
		line2_no_comma = sprintf('$s=%i$, $h=%0.3f$, $N=%i$', s, h, N);

		title_str_iter_part = 'Iterations Required at Each Step';
		title_str_iter = sprintf('%s\n%s', title_str_iter_part, method_str, line2);

		figure
		plot(1:N, iter_vec, 'bo', 'LineWidth', 1.5)
		hold on
		xlabel('Step', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('Iterations', 'Interpreter', 'latex', 'FontSize', 16)
		set(gca,'FontSize',16)
		title(title_str_iter,'Interpreter','latex','FontSize',24)
		grid on
		hold off
	end

	% --------------------
	% Other plots
	% --------------------
	if strcmpi(plot_str, 'pendulum')

		line2 = sprintf(', $s=%i$, $h=%0.3f$, $N=%i$', s, h, N);

		% Plot the position coordinates
		titlestr0 = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
		g_title_str = sprintf('%s\n%s', titlestr0, method_str, line2);
		t = linspace(0,2*pi);

		figure
		plot(yns(1,:),yns(2,:), 'bo', 'LineWidth', 1.5)
		hold on
		plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
		xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(g_title_str,'Interpreter','latex','FontSize',24)
		xlim([-1 1])
		ylim([-1 1])
		axis square
		grid on
		legend('Numerical Solution', 'The Manifold $y^Ty-l^2=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
		hold off

		% plot 2y^Tz vs. time
		for i=1:length(yns)
			if i==1
				values = 2*dot(yns(:,i), zns(:,i));
			else
				values = [values 2*dot(yns(:,i), zns(:,i))];
			end
		end

		line1 = 'Adherence of the Solution to $2 y^T z=0$';

		gyf_title_str = sprintf('%s\n%s%s', line1, method_str, line2);

		figure
		plot(0:h:N*h, values, 'b-', 'LineWidth', 1.5)
		hold on
		%plot(0:h:N*h, zeros(1,N+1), 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$2y(t)^Tz(t)$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(gyf_title_str,'Interpreter','latex','FontSize',24)
		ylim([-0.2 0.2])
		axis square
		grid on
		legend('Numerical Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

		% Plot y1 vs t, y2 vs t, and energy vs t

		titlestr1 = '$y_1$ vs.~$t$';
		titlestr2 = '$y_2$ vs.~$t$';
		titlestr3 = 'Energy vs.~$t$';

		figure
		tiledlayout(3,1)
		nexttile
		plot(0:h:N*h, yns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr1,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, yns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr2,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr3,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		% Plot energy vs. t in a separate figure
		figure
		plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr3,'Interpreter','latex','FontSize',24)
		grid on
		hold off

	elseif strcmpi(plot_str, 'pendulum theta')

		% change notation to (theta,omega)
		thetans = yns;
		omegans = zns;

		% find position and velocity coordinates (q1,q2) and (p1,p2)
		qns(1,:) = l*sin(thetans);
		qns(2,:) = -l*cos(thetans);

		% note this uses component-wise multiplication
		pns(1,:) = l*omegans.*cos(thetans); 
		pns(2,:) = l*omegans.*sin(thetans); 
		
		% plot theta(t) and omega(t)
		figure
		plot(0:h:N*h, thetans, 'b-', 'LineWidth', 1.5)
		hold on
		plot(0:h:N*h, omegans, 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		%ylabel(labels_y{j}, 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title('$\theta$ and $\omega$ vs.~$t$','Interpreter','latex','FontSize',24)
		legend('$\theta(t)$', '$\omega(t)$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		grid on
		hold off
		
		% plot q1, q2, p1, p2
		figure
		tiledlayout(4,1)
		nexttile
		plot(0:h:N*h, qns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		title('Position and Velocity Coordinates', 'Interpreter', 'latex', 'FontSize', 24)
		hold off

		nexttile
		plot(0:h:N*h, qns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, pns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$z_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, pns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$z_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		% plot position coordinates on constraint manifold
		line2 = sprintf(', $s=%i$, $h=%0.3f$, $N=%i$', s, h, N);

		% Plot the position coordinates
		titlestr0 = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
		g_title_str = sprintf('%s\n%s', titlestr0, method_str, line2);
		t = linspace(0,2*pi);

		figure
		plot(qns(1,:),qns(2,:), 'bo', 'LineWidth', 1.5)
		hold on
		plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
		xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(g_title_str,'Interpreter','latex','FontSize',24)
		xlim([-1 1])
		ylim([-1 1])
		axis square
		grid on
		legend('Numerical Solution', 'The Manifold $y^Ty-l^2=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
		hold off

		% plot 2y^Tz vs. time
		for i=1:length(yns)
			if i==1
				values = 2*dot(qns(:,i), pns(:,i));
			else
				values = [values 2*dot(qns(:,i), pns(:,i))];
			end
		end

		line1 = 'Adherence of the Solution to $2 y^T z=0$';

		gyf_title_str = sprintf('%s\n%s%s', line1, method_str, line2);

		figure
		plot(0:h:N*h, values, 'b-', 'LineWidth', 1.5)
		hold on
		%plot(0:h:N*h, zeros(1,N+1), 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$2y(t)^Tz(t)$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(gyf_title_str,'Interpreter','latex','FontSize',24)
		ylim([-0.2 0.2])
		axis square
		grid on
		legend('Numerical Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off
		
	elseif strcmpi(plot_str, 'double pendulum angles')

		% change notation to (theta,omega)
		thetans = yns;
		omegans = zns;

		% find position and velocity coordinates (q1,q2,q3,q4) and (p1,p2,p3,p4)
		qns(1,:) = l1*sin(thetans(1,:));
		qns(2,:) = -l1*cos(thetans(1,:));
		qns(3,:) = l1*sin(thetans(1,:)) + l2*sin(thetans(2,:));
		qns(4,:) = -l1*cos(thetans(1,:)) - l2*cos(thetans(2,:));

		% note this uses component-wise multiplication
		pns(1,:) = l1*omegans(1,:).*cos(thetans(1,:)); 
		pns(2,:) = l1*omegans(1,:).*sin(thetans(1,:)); 
		pns(3,:) = l1*omegans(1,:).*cos(thetans(1,:)) + l2*omegans(2,:).*cos(thetans(2,:)); 
		pns(4,:) = l1*omegans(1,:).*sin(thetans(1,:)) + l2*omegans(2,:).*sin(thetans(2,:)); 
		
		% plot theta(t) and omega(t)
		figure
		plot(0:h:N*h, thetans(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		plot(0:h:N*h, thetans(2,:), 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$\theta_i$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title('$\theta$ vs.~$t$','Interpreter','latex','FontSize',24)
		legend('$\theta_1(t)$', '$\theta_2(t)$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		grid on
		hold off

		figure
		plot(0:h:N*h, omegans(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		plot(0:h:N*h, omegans(2,:), 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$\omega_i$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title('$\omega$ vs.~$t$','Interpreter','latex','FontSize',24)
		legend('$\omega_1(t)$', '$\omega_2(t)$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		grid on
		hold off
		
		% plot q1, q2, q3, q4
		figure
		tiledlayout(4,1)
		nexttile
		plot(0:h:N*h, qns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		title('Position Coordinates', 'Interpreter', 'latex', 'FontSize', 24)
		hold off

		nexttile
		plot(0:h:N*h, qns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, qns(3,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_3$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, qns(4,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_4$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		% plot p1, p2, p3, p4
		figure
		tiledlayout(4,1)
		nexttile
		plot(0:h:N*h, pns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		title('Velocity Coordinates', 'Interpreter', 'latex', 'FontSize', 24)
		hold off

		nexttile
		plot(0:h:N*h, pns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, pns(3,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_3$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, pns(4,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_4$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		%title(title_str_y{j},'Interpreter','latex','FontSize',24)
		grid on
		hold off

		% plot deviation from constraints
		c1_vals = qns(1,:).^2 + qns(2,:).^2 - l1^2;
		c2_vals = (qns(1,:)-qns(3,:)).^2 + (qns(2,:)-qns(4,:)).^2 - l2^2;
		c3_vals = qns(1,:).*pns(1,:) + qns(2,:).*pns(2,:);
		c4_vals = (qns(1,:)-qns(3,:)).*(pns(1,:)-pns(3,:)) + (qns(2,:)-qns(4,:)).*(pns(2,:)-pns(4,:));

		line1 = 'Adherence of the Solution to Constraints';
		constraint_title_str = sprintf('%s\n%s%s', line1, method_str, line2);

		figure
		plot(0:h:N*h, c1_vals, 'b-', 'LineWidth', 1.5)
		hold on
		plot(0:h:N*h, c2_vals, 'r-', 'LineWidth', 1.5)
		plot(0:h:N*h, c3_vals, 'k-', 'LineWidth', 1.5)
		plot(0:h:N*h, c4_vals, 'c-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$c_i(q,p)$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(constraint_title_str,'Interpreter','latex','FontSize',24)
		%ylim([-0.2 0.2])
		axis square
		grid on
		legend('$c_1(q(t),p(t))$', '$c_2(q(t),p(t))$', '$c_3(q(t),p(t))$', '$c_4(q(t),p(t))$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

		% plot position coordinates
		figure
		plot(qns(1,:), qns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		plot(qns(3,:), qns(4,:), 'r-', 'LineWidth', 1.5)
		% initial points
		h3 = plot(qns(1,1), qns(2,1), 'bo', 'MarkerFaceColor','b');
		h4 = plot(qns(3,1), qns(4,1), 'ro', 'MarkerFaceColor','r');
		xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title('Cartesian Coordinates of the Double Pendulum','Interpreter','latex','FontSize',24)
		%ylim([-0.2 0.2])
		axis square
		grid on
		legend('Mass 1', 'Mass 2', 'Mass 1 Initial', 'Mass 2 Initial', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEastoutside')
		hold off

	elseif strcmpi(plot_str, 'linear')

		% plot deviation from constraints
		c_vals = yns(1,:).*yns(2,:).^2-1;

		line1 = sprintf('Adherence of the Solution to Constraint\n(Partitioned ODE), ');
		constraint_title_str = sprintf('%s%s\n%s', line1, method_str, line2_no_comma);

		figure
		plot(0:h:N*h, c_vals, 'r-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$c(t)$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(constraint_title_str,'Interpreter','latex','FontSize',24)
		ax=gca;
		ax.YAxis.Exponent = 0;     % turn off the top "×10^n"
		ytickformat('%.1e')   % if there's an exponent, it goes at every tick mark
		%ylim([-0.2 0.2])
		axis square
		grid on
		%legend('$c_1(q(t),p(t))$', '$c_2(q(t),p(t))$', '$c_3(q(t),p(t))$', '$c_4(q(t),p(t))$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

	elseif strcmpi(plot_str, 'two body')

		line2 = sprintf('%s, $s=%i$, $h=%0.3f$, $N=%i$', method_str, s, h, N);
		
		% Plot the position coordinates
		figure
		plot(yns(1,:),yns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		
		line0 = sprintf('Position Coordinates $(y_1,y_2)$');
		line1 = sprintf('of the Two-Body Problem');
		title_str = sprintf('%s\n%s\n%s', line0, line1, line2);

		% plot the analytical solution
		theta = linspace(0,2*pi);
		radius = alpha0 ./ (1 + e*cos(theta));
		x1_true = radius .* cos(theta);
		x2_true = radius .* sin(theta);
		plot(x1_true,x2_true,'r--','LineWidth',2)

		xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(title_str,'Interpreter','latex','FontSize',24)
		xlim([-2.2 1.2])
		ylim([-1.5 2])
		xticks([-2 -1 0 1])
		yticks([-1 0 1])
		axis square
		grid on
		legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

		% Make the energy vector.
		% Note that here energy E = H(q,p).
		qns = yns;
		pns = zns;

		for i=1:N+1
			qn = qns(:,i);
			pn = pns(:,i);
			if i==1
				E = (1/2)*norm(pn)^2 - 1/norm(qn);
			else
				Enp = (1/2)*norm(pn)^2 - 1/norm(qn);
				E = [E Enp];
			end
		end

		titlestr3 = 'Energy vs.~$t$';
		figure
		plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr3,'Interpreter','latex','FontSize',24)
		grid on
		hold off

	elseif strcmpi(plot_str, 'lotka-volterra')

		line2 = sprintf('%s, $s=%i$, $h=%0.3f$, $N=%i$', method_str, s, h, N);
		
		% Plot the phase coordinates
		figure
		plot(yns(1,:),zns(1,:), 'bo-', 'LineWidth', 1.5)
		hold on
		
		title_str = sprintf('%s\n%s', 'The Lotka-Volterra Model', line2);

		% plot the analytical solution
		% Note this assumes the constants c0=2, a0=1, gamma0=1, alpha0=1

		I0 = 3*log(2)-4;
		I = @(y,z) log(y) - y + 2*log(z) - z - I0;
		fimplicit(I, [0.01 5 0.01 5], 'r--','LineWidth',2)

		xlabel('$y$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$z$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(title_str,'Interpreter','latex','FontSize',24)
		xlim([0 4])
		ylim([0 4])
		%xticks([-2 -1 0 1])
		%yticks([-1 0 1])
		axis square
		grid on
		legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

	end

	% ------------------------------------------------------------------------------------------
	% Plot components for all plot strings except 'pendulum theta' and 'double pendulum angles'
	% ------------------------------------------------------------------------------------------

	if ~strcmpi(plot_str, 'pendulum theta') && ~strcmpi(plot_str, 'double pendulum angles')
		% Plot y_j vs t and z_j vs t for all j.

		% If ny+nz is below or equal to this threshold, all plots appear in
		% a single figure.  Past this threshold, the y, z, and u components
		% are plotted separately.
		tiled_layout_threshold = 5;

		title_str_y = {};
		title_str_z = {};

		labels_y = {};
		labels_z = {};

		for j=1:ny
			title_str_y{j} = sprintf('$y_%i$ vs.~$t$',j);
			labels_y{j} = sprintf('$y_%i$',j);
		end

		for j=1:nz
			title_str_z{j} = sprintf('$z_%i$ vs.~$t$',j);
			labels_z{j} = sprintf('$z_%i$',j);
		end

		figure

		if ny+nz <= tiled_layout_threshold
			tiledlayout(ny+nz,1)
		else
			tiledlayout(ny,1)
		end

		for j=1:ny
			nexttile
			plot(0:h:N*h, yns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_y{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_y{j},'Interpreter','latex','FontSize',24)
			grid on
			hold off
		end

		if ny+nz > tiled_layout_threshold 
			figure
			tiledlayout(nz,1)
		end

		for j=1:nz
			nexttile
			plot(0:h:N*h, zns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_z{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_z{j},'Interpreter','latex','FontSize',24)
			grid on
			hold off
		end
	end
end

