clear;

T = 1:1000;
lambda = repmat([50], 1, length(T));
lambda(length(T)/2) = 60;
mu1 =  horzcat((repmat([25], 1, length(T)/2)), (repmat([35], 1, length(T)/2)));
mu2 =  horzcat((repmat([25], 1, length(T)/2)), (repmat([15], 1, length(T)/2)));

mu = vertcat(mu1, mu2);

function qdot = q(t, x, lamb, mu)

	den = (1 + x(1)) * mu(1) + (1 + x(2)) * mu(2);
	f1 = lamb * (1 + x(1)) * mu(1) / den;
	f2 = lamb * (1 + x(2)) * mu(2) / den;

	qd1 = f1 - mu(1);
	if (qd1 < 0)
		qd1 = f1;
	endif

	qd2 = f2 - mu(2);
	if (qd2 < 0)
		qd2 = f2;
	endif

	qdot(1) = f1 - mu(1);
	qdot(2) = f2 - mu(2);

endfunction

x0_1 = 200;
x0_2 = 200;

res_start_x = [];
res_start_y = [];

for i = 1:(length(T)-1)
	res = ode23(@q, T(i:i+1), [x0_1; x0_2], lambda(i), [mu(1, i); mu(2, i)]);
	x0_1 = res.y(length(res.y), 1);
	x0_2 = res.y(length(res.y), 2);
	if (i == 1)
		res_start_x = res.x;
		res_start_y = res.y;
	else
		res_start_x = vertcat(res_start_x, res.x);
		res_start_y = vertcat(res_start_y, res.y);
	endif
endfor


plot(res_start_x, res_start_y)