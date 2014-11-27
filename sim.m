clear;
pkg load odepkg;

Tstart = 1;
Tend = 600;
NumT = Tend * 2;
T = linspace(Tstart, Tend, NumT);

function qdot = q(t, x, xd,  num_clients, num_servers, lags, qsz_exponent, os)

	arrival_rate = ArrivalRate(t, num_clients);
	flow_matrix = zeros(num_clients, num_servers);
	for client = 1:num_clients
		qmu_inverse = zeros(num_servers, 1);
		for server = 1:num_servers
			measured_q = max(xd(server, server), 0);
			% 1/rate == mu
			% need a way to import time lag
			os_n = os * num_clients;
			service_rate = ServiceRate(t - lags(server), num_servers)(server);
			qmu_inverse(server) = service_rate/((1 + os_n + measured_q) ^ qsz_exponent);
		endfor
		TotalWeight = sum(qmu_inverse);

		% xxx: use vector operation
		for server = 1:num_servers
			flow_fraction = qmu_inverse(server)/(TotalWeight);
			flow_matrix(client, server) = arrival_rate(client) * flow_fraction;
		endfor
	endfor
	flow_matrix;

	% use the current rate to find the integral
	rate = ServiceRate(t, num_servers);

	% printf("%d %d %d %d %d %d %d %d\n", x(1), x(2), t, f1, f2, xd(2, 2), xd(1, 1), lamb);

	for i = 1:num_servers
		if (x(i) <= 0)
			rate(i) = 0;
		endif
		qdot(i) = sum(flow_matrix(:, i)) - rate(i);
	endfor
endfunction

function ArrivalRate = ArrivalRate(t, num_clients)
	lamb = 30;
	ArrivalRate = repmat(lamb, 1, num_clients);
	% ArrivalRate = [20; 20; 20];
endfunction

% xxx: use sinusoid
function ServiceRate = ServiceRate(t, num_servers)
	
	ServiceRate = repmat(50, 1, num_servers);

	if(t > 400)
		ServiceRate = [90;50;10];
		% ServiceRate = [30;70;50];
	elseif(t > 100)
		ServiceRate = [30;70;50];
	endif

	% st = 50;
	% amplitude = 5;
	% ServiceRate = [sin(t/10 + pi/4) * amplitude + st; 
	% 			   sin(t/10 + pi/2) * amplitude + st; 
	% 			   sin(t/10 + 3* pi/4) * amplitude + st];
endfunction

initq = 500;
num_clients = 5;
num_servers = 3;
qsz_exponent = 1;
os = 0;

init = repmat([initq], num_servers, 1);
lags = repmat([0.5], 1, num_servers);

hist_mat = repmat(initq, num_servers, num_servers);

res = ode23d(@q, T, init, lags, hist_mat, num_clients, num_servers, lags, qsz_exponent, os);

% XXX: Add legend
plot(res.x, res.y, 'LineWidth', 2);
