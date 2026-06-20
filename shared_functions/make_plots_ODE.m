function make_plots_ODE(yns,iter_vec,plots_params)
	% function make_plots_ODE(yns,iter_vec,plots_params)
	%
	% Make plots for ODEs.
	%
	% Inputs:
	%	yns: The matrices consisting of the update steps yn
	%	iter_vec: vector of iterations 
	%	plots_params: Structured array with the following fields:
	%		h: stepsize
	%		N: number of steps 
	%		s: number of stages
	%		problem_str: describes the DAE to be solved.
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

	if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1') || strcmpi(problem_str, 'pendulum version 2')
		E = plots_params.E; % energy
	elseif strcmpi(problem_str, 'two body')
		alpha0 = plots_params.two_body_struct.alpha0;
		T = plots_params.two_body_struct.T;
		a = plots_params.two_body_struct.a;
		e = plots_params.two_body_struct.e;
	end

	ny = size(yns,1);

	% --------------------
	% Plots of iterations 
	% --------------------
	if plot_iterations

		line2 = sprintf(', $s=%i$, $h=%0.3f$, $N=%i$', s, h, N);

		title_str_iter_part = 'Iterations Required at Each Step';
		title_str_iter = sprintf('%s\n%s', title_str_iter_part, method_str, line2);

		figure
		plot(1:N, iter_vec, 'bo-', 'LineWidth', 1.5)
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
		plot(yns(1,:),yns(2,:), 'bo-', 'LineWidth', 1.5)
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
		legend('Numerical Solution', 'The Manifold $y_1^2+y_2^2=l^2$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
		hold off

		% plot 2y^Tz vs. time
		for i=1:length(yns)
			if i==1
				values = 2*dot(yns(1:2,i), yns(3:4,i));
			else
				values = [values 2*dot(yns(1:2,i), yns(3:4,i))];
			end
		end

		% for partitioned ODEs and DAEs, this would read '2y^Tz=0'
		line1 = 'Adherence of the Solution to $2 (y_1 y_3 + y_2 y_4)=0$';

		gyf_title_str = sprintf('%s\n%s%s', line1, method_str, line2);

		figure
		plot(0:h:N*h, values, 'b-', 'LineWidth', 1.5)
		hold on
		plot(0:h:N*h, zeros(1,N+1), 'r-', 'LineWidth', 1.5)
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$2(y_1y_3 + y_2y_4)$', 'Interpreter', 'latex', 'FontSize', 16)
		%set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(gyf_title_str,'Interpreter','latex','FontSize',24)
		ylim([-0.2 0.2])
		axis square
		grid on
		legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
		hold off

		% Plot energy vs. t in a separate figure
		figure
		plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title('Energy vs.~$t$','Interpreter','latex','FontSize',24)
		grid on
		hold off

	elseif strcmpi(plot_str, 'two body')

		line2 = sprintf('%s, $s=%i$, $h=%0.3f$, $N=%i$', method_str, s, h, N);
		
		% Plot the position coordinates
		figure
		plot(yns(1,:),yns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		
		title_str = sprintf('%s\n%s', 'Pos.~Coordinates $(y_1,y_2)$ of the Two-Body Problem', line2);

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
		qns = yns(1:2,:);
		pns = yns(3:4,:);

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
	end

	% ------------------------------------------------------
	% Plot components, no matter what the plot_str value is.
	% ------------------------------------------------------

	% Plot y_j vs t for all j.
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

	figure

	tiledlayout(ny,1)

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
end

