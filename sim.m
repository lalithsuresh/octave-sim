clear;
pkg load odepkg;

T = 1:1000;

function qdot = q(t, x, xd)
	lamb = lambda(t);
	mu = service_rate(t);
	measured_q1 = max(xd(1, 1), 0);
	measured_q2 = max(xd(2, 2), 0);

	qmu_1 = (1 + measured_q1)/mu(1);
	qmu_2 = (1 + measured_q2)/mu(2);
	C = qmu_2/qmu_1;

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

	printf("%d %d %d %d %d %d %d %d\n", x(1), x(2), t, f1, f2, xd(2, 2), xd(1, 1), lamb);
	qdot(1) = f1 - mu_1;
	qdot(2) = f2 - mu_2;
endfunction

function lambda = lambda(t)
	lamb = 50;

	if (t > 30)
		lamb = 50;
	elseif (t > 20)
		lamb = 20;
	endif

	lambda(1) = lamb;
endfunction

function service_rate = service_rate(t)

	st_1 = 35;
	st_2 = 15;

	if (t > 30)
		st_1 = 15;
		st_2 = 35;
	elseif (t > 20)
		st_1 = 35;
		st_2 = 25;
	endif

	service_rate(1) = st_1;
	service_rate(2) = st_2;

endfunction

initq = 100;
x0_1 = initq;
x0_2 = initq;

hist_mat = repmat(initq, 2, 3);

res = ode23d(@q, T, [x0_1; x0_2], [1, 1], hist_mat);

plot(res.x, res.y)
