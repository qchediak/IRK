function make_plots(yns,zns,uns,iterations_struct,plots_params)
	% function make_plots(yns,zns,uns,iterations_struct,plots_params)
	%
	% Make plots for index 3 DAEs.
	%
	% Inputs:
	%	yns, zns, uns: The matrices consisting of the update steps (yn,zn,un)
	%	iterations_struct: Structured array of vectors of iterations 
	%	plots_params: Structured array with the following fields:
	%		h: stepsize
	%		N: number of steps 
	%		s: number of stages
	%		problem_str: describes the DAE to be solved.
	%		method_str: describes the PRK method
	%		plot_str: sets a plotting option.
	%		iterative_scheme_struct: structured array of strings
	%			describing the iterative method(s)
	%		number_nonlinear_systems: The number of nonlinear systems. 
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
	iterative_scheme_struct = plots_params.iterative_scheme_struct;
	number_nonlinear_systems = plots_params.number_nonlinear_systems;
	plot_iterations = plots_params.plot_iterations; 

	% unpack the true solutions
	if strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'linear') || strcmpi(problem_str, 'second nonlinear')
		y_true = plots_params.true_solution.y_true;
		z_true = plots_params.true_solution.z_true;
		u_true = plots_params.true_solution.u_true;
	elseif strcmpi(problem_str, 'pendulum')
		E = plots_params.E; % energy
	end

	ny = size(yns,1);
	nz = size(zns,1);
	nu = size(uns,1);

	% --------------------
	% Plots of iterations 
	% --------------------
	if plot_iterations
		% unpack the structured array iterations_struct 
		if number_nonlinear_systems==1
			iter_vec = iterations_struct.iter_vec;
		elseif number_nonlinear_systems==2
			iter1_vec = iterations_struct.iter1_vec;
			iter2_vec = iterations_struct.iter2_vec;
		end

		if number_nonlinear_systems==1

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

		elseif number_nonlinear_systems==2

			line2 = sprintf(', $s=%i$, $h=%0.3f$, $N=%i$', s, h, N);

			title_str_iter1_part = 'Iterations Required at Each Step (Nonlinear System 1)';
			title_str_iter1 = sprintf('%s\n%s', title_str_iter1_part, method_str, line2);

			title_str_iter2_part = 'Iterations Required at Each Step (Nonlinear System 2)';
			title_str_iter2 = sprintf('%s\n%s', title_str_iter2_part, method_str, line2);

			figure
			tiledlayout(2,1)
			nexttile
			plot(1:N, iter1_vec, 'bo', 'LineWidth', 1.5)
			hold on
			xlabel('Step', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel('Iterations', 'Interpreter', 'latex', 'FontSize', 16)
			set(gca,'FontSize',16)
			title(title_str_iter1,'Interpreter','latex','FontSize',24)
			grid on
			hold off

			nexttile
			plot(1:N, iter2_vec, 'bo', 'LineWidth', 1.5)
			hold on
			xlabel('Step', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel('Iterations', 'Interpreter', 'latex', 'FontSize', 16)
			set(gca,'FontSize',16)
			title(title_str_iter2,'Interpreter','latex','FontSize',24)
			grid on
			hold off

		else

			error('Error (make_plots): Invalid value of number_nonlinear_systems')

		end
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
		legend('Numerical Solution', 'The Manifold $g(y)=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
		hold off

		% plot 2y^Tz vs. time
		for i=1:length(yns)
			if i==1
				values = 2*dot(yns(:,i), zns(:,i));
			else
				values = [values 2*dot(yns(:,i), zns(:,i))];
			end
		end

		line1 = 'Adherence of the Solution to $(g_y f)(y,z)=0$';

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

	elseif strcmpi(plot_str, 'linear')

		% Plot constraint violation
		% plot deviation from constraints
		c_vals = yns(1,:).*yns(2,:).^2-1;

		line1 = sprintf('Adherence of the Solution to Constraint (DAE)');
		constraint_title_str = sprintf('%s\n%s%s', line1, method_str, line2);

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

		% plot components and the true solution
		title_str_y = {};
		title_str_z = {};
		title_str_u = {};

		labels_y = {};
		labels_z = {};
		labels_u = {};

		for j=1:ny
			title_str_y{j} = sprintf('$y_%i$ vs.~$t$',j);
			labels_y{j} = sprintf('$y_%i$',j);
		end

		for j=1:nz
			title_str_z{j} = sprintf('$z_%i$ vs.~$t$',j);
			labels_z{j} = sprintf('$z_%i$',j);
		end

		for j=1:nu
			title_str_u{j} = sprintf('$u_%i$ vs.~$t$',j);
			labels_u{j} = sprintf('$u_%i$',j);
		end

		figure

		tiledlayout(ny+nz+nu,1)

		% values of the true solution
		y_true_vals = y_true(0:h:N*h);
		z_true_vals = z_true(0:h:N*h);
		u_true_vals = u_true(0:h:N*h);

		for j=1:ny
			nexttile
			plot(0:h:N*h, yns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, y_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_y{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_y{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end
		
		for j=1:nz
			nexttile
			plot(0:h:N*h, zns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, z_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_z{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_z{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end

		for j=1:nu
			nexttile
			plot(0:h:N*h, uns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, u_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_u{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_u{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end

	elseif strcmpi(plot_str, 'component') && (strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'linear') || strcmpi(problem_str, 'second nonlinear'))
		% For these problems, also plot the true solution
		title_str_y = {};
		title_str_z = {};
		title_str_u = {};

		labels_y = {};
		labels_z = {};
		labels_u = {};

		for j=1:ny
			title_str_y{j} = sprintf('$y_%i$ vs.~$t$',j);
			labels_y{j} = sprintf('$y_%i$',j);
		end

		for j=1:nz
			title_str_z{j} = sprintf('$z_%i$ vs.~$t$',j);
			labels_z{j} = sprintf('$z_%i$',j);
		end

		for j=1:nu
			title_str_u{j} = sprintf('$u_%i$ vs.~$t$',j);
			labels_u{j} = sprintf('$u_%i$',j);
		end

		figure

		tiledlayout(ny+nz+nu,1)

		% values of the true solution
		y_true_vals = y_true(0:h:N*h);
		z_true_vals = z_true(0:h:N*h);
		u_true_vals = u_true(0:h:N*h);

		for j=1:ny
			nexttile
			plot(0:h:N*h, yns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, y_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_y{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_y{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end
		
		for j=1:nz
			nexttile
			plot(0:h:N*h, zns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, z_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_z{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_z{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end

		for j=1:nu
			nexttile
			plot(0:h:N*h, uns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			plot(0:h:N*h, u_true_vals(j,:), 'r--', 'LineWidth', 1.5)
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_u{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_u{j},'Interpreter','latex','FontSize',24)
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'eastoutside')
			hold off
		end

	elseif strcmpi(plot_str, 'component') 

		% Plot y_j vs t, z_j vs t, and u_j vs t for all j.

		% If ny+nz+nu is below or equal to this threshold, all plots appear in
		% a single figure.  Past this threshold, the y, z, and u components
		% are plotted separately.
		tiled_layout_threshold = 5;

		title_str_y = {};
		title_str_z = {};
		title_str_u = {};

		labels_y = {};
		labels_z = {};
		labels_u = {};

		for j=1:ny
			title_str_y{j} = sprintf('$y_%i$ vs.~$t$',j);
			labels_y{j} = sprintf('$y_%i$',j);
		end

		for j=1:nz
			title_str_z{j} = sprintf('$z_%i$ vs.~$t$',j);
			labels_z{j} = sprintf('$z_%i$',j);
		end

		for j=1:nu
			title_str_u{j} = sprintf('$u_%i$ vs.~$t$',j);
			labels_u{j} = sprintf('$u_%i$',j);
		end

		figure

		if ny+nz+nu <= tiled_layout_threshold
			tiledlayout(ny+nz+nu,1)
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

		if ny+nz+nu > tiled_layout_threshold 
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

		if ny+nz+nu > tiled_layout_threshold
			figure
			tiledlayout(nu,1)
		end

		for j=1:nu
			nexttile
			plot(0:h:N*h, uns(j,:), 'b-', 'LineWidth', 1.5)
			hold on
			xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel(labels_u{j}, 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(title_str_u{j},'Interpreter','latex','FontSize',24)
			grid on
			hold off
		end

	else

		error('Error (make_plots): Invalid plot_str and problem_str combination')

	end
end

