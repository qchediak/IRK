function make_movie_ODE(yns,plots_params)
	% function make_movie_ODE(yns,plots_params)
	%
	% Make a movie and play it. Specific to the pendulum problem.
	%
	% Inputs:
	%	yns: Matrix of the update steps
	%	plots_params: Structured array of parameters.
	%
	% Outputs:
	%	None.

	h = plots_params.h;
	N = plots_params.N;
	problem_str = plots_params.problem_str;

	if strcmpi(problem_str, 'two body')
		alpha0 = plots_params.two_body_struct.alpha0;
		e = plots_params.two_body_struct.e;
	end	

	if strcmpi(problem_str, 'pendulum')
		titlestr0 = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
		t = linspace(0,2*pi);

		moviein(N);
		%figure
		for j=1:length(yns)
			plot(yns(1,j),yns(2,j), 'bo', 'LineWidth', 3.0)
			hold on
			plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
			xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(titlestr0,'Interpreter','latex','FontSize',24)
			xlim([-1 1])
			ylim([-1 1])
			axis square
			grid on
			legend('Numerical Solution', 'The Manifold $g(y)=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
			hold off
			
			M(:,j) = getframe;
		end

	elseif strcmpi(problem_str, 'two body')

		titlestr0 = 'Pos.~Coordinates $(y_1,y_2)$ of the Two-Body Problem';
		t = linspace(0,2*pi);

		theta = linspace(0,2*pi);
		radius = alpha0 ./ (1 + e*cos(theta));
		x1_true = radius .* cos(theta);
		x2_true = radius .* sin(theta);

		moviein(N);
		%figure
		for j=1:length(yns)
			plot(yns(1,j),yns(2,j), 'bo', 'LineWidth', 3.0)
			hold on
			plot(x1_true,x2_true,'r--','LineWidth',2)
			xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
			ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
			set(get(gca,'YLabel'),'Rotation',0)
			set(gca,'FontSize',16)
			title(titlestr0,'Interpreter','latex','FontSize',24)
			xlim([-2.2 1.2])
			ylim([-1.5 2])
			xticks([-2 -1 0 1])
			yticks([-1 0 1])
			axis square
			grid on
			legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
			hold off
			
			M(:,j) = getframe;
		end
	end

	% play the movie
	movie(M,1,10)
end


