clear;
pkg load odepkg;

T = 1:100;
mu1 =  horzcat((repmat([25], 1, length(T)/2)), (repmat([45], 1, length(T)/2)));
mu2 =  horzcat((repmat([25], 1, length(T)/2)), (repmat([15], 1, length(T)/2)));

mu = vertcat(mu1, mu2);

function qdot = q(t, x, xd)
	lamb = lambda(t);
	mu = service_time(t);
	qmu_1 = (1 + xd(1, 1)) * mu(1);
	qmu_2 = (1 + xd(2, 2)) * mu(2);
	C = qmu_1/qmu_2;

	% den = qmu_1 + qmu_2;
	% f1 = lamb * qmu_1 / den;
	% f2 = lamb * qmu_2 / den;


	% C = (mu(1)/mu(2)) ^ (1/3);
	% f1 = lamb/C;
	% f2 = lamb - f1;
	% f1 = (C * lamb + C * xd(2, 2) - xd(1, 1))/(C + 1);
	% f2 = lamb - f1;
	f1 = C * lamb / (1 + C);
	f2 = lamb - f1;

	mu_1 = mu(1);
	if (x(1) <= 0)
		mu_1 = 0;
	endif

	mu_2 = mu(2);
	if (x(2) <= 0)
		mu_2 = 0;
	endif

	% printf("%d %d %d %d %d %d %d\n", x(1), x(2), t, f1, f2, xd(2, 2), xd(1, 1));
	qdot(1) = f1 - mu_1;
	qdot(2) = f2 - mu_2;
endfunction

function lambda = lambda(t)
	lambda(1) = 50;
endfunction

function service_time = service_time(t)
	service_time(1) = 25;
	service_time(2) = 25;
endfunction

initq = 1;
x0_1 = initq;
x0_2 = initq;

res_start_x = [];
res_start_y = [];

hist_mat = repmat(initq, 2, 3);

% for i = 1:1:(length(T)-1)
res = ode23d(@q, T, [x0_1; x0_2], [0.01, 0.01], hist_mat);
% 	x0_1 = res.y(length(res.y), 1);
% 	x0_2 = res.y(length(res.y), 2);
% 	if (i == 1)
% 		res_start_x = res.x;
% 		res_start_y = res.y;
% 	else
% 		res_start_x = vertcat(res_start_x, res.x);
% 		res_start_y = vertcat(res_start_y, res.y);
% 	endif
% 	hist_mat = (res.y((length(res.y) - 1):length(res.y),:)).';
% 	% hist_mat = res
% endfor

plot(res.x, res.y)
% plot(res_start_x, res_start_y)
