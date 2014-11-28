clear;
pkg load odepkg;

Tstart = 1;
Tend = 600;
NumT = Tend;
T = linspace(Tstart, Tend, NumT);
% T = [0,600];

function qdot = q(t, x, xd,  num_clients, num_servers, lags, qsz_exponent, os)
	qdot = zeros(num_servers + num_clients, 1);
	arrival_rate = ArrivalRate(t, num_clients);
	flow_matrix = zeros(num_clients, num_servers);
	sending_rate = [10; 10; 10];
	% if(t > 402)
	% 	sending_rate = [56/num_clients; 10; 10];
	% endif
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
		backlogRate = 0;
		total_to_allocate = arrival_rate(client);
		is_full = zeros(num_servers, 1);

		if (x(num_servers + client) > 0)
			% When in backpressure mode, fire away at full rate
			total_to_allocate = sum(sending_rate);
		endif

		for k = 1:num_servers
			for server = 1:num_servers
				if (is_full(server) == 0)
					flow_fraction = qmu_inverse(server)/(TotalWeight);
					flow_matrix(client, server) = total_to_allocate * flow_fraction;

					if(flow_matrix(client, server) > sending_rate)
						backlogRate = backlogRate + flow_matrix(client, server) - sending_rate(server);
						flow_matrix(client, server) = sending_rate(server);
					endif

					if (flow_matrix(client, server) >= sending_rate)
						is_full(server) = 1;
						TotalWeight = TotalWeight - qmu_inverse(server);
					endif
				endif
			endfor

			if (backlogRate > 0)
				total_to_allocate = backlogRate;
			else
				break
			endif
		endfor

		% Update backlog size at client
		total_drain = sum(flow_matrix(client, :));
		qdot(num_servers + client) = arrival_rate(client) - total_drain;
	endfor
	% if(t > 400)
	% 	printf("LOL %d\n", t);
	% 	flow_matrix
	% endif

	% use the current rate to find the integral
	rate = ServiceRate(t, num_servers);

	% printf("%d %d %d %d %d %d %d %d\n", x(1), x(2), t, f1, f2, xd(2, 2), xd(1, 1), lamb);

	for server = 1:num_servers
		if (x(server) <= 0)
			rate(server) = 0;
		endif
		qdot(server) = sum(flow_matrix(:, server)) - rate(server);
	endfor
endfunction

function ArrivalRate = ArrivalRate(t, num_clients)
	lamb = 30;
	ArrivalRate = repmat(lamb, 1, num_clients);

	% if(t > 401)
	% 	ArrivalRate = repmat(lamb, 1, num_clients);
	% elseif(t > 400)
	% 	ArrivalRate = repmat(31, 1, num_clients);
	% endif

endfunction

% xxx: use sinusoid
function ServiceRate = ServiceRate(t, num_servers)
	
	ServiceRate = repmat(50, 1, num_servers);

	% if(t > 401)
	% 	ServiceRate = [56;50;50];
	% 	% ServiceRate = [30;70;50];
	% elseif(t > 100)
	% 	% ServiceRate = [30;70;50];
	% 	ServiceRate = [50;50;50];
	% endif

	% st = 50;
	% amplitude = 5;
	% ServiceRate = [sin(t/10 + pi/4) * amplitude + st; 
	% 			   sin(t/10 + pi/2) * amplitude + st; 
	% 			   sin(t/10 + 3* pi/4) * amplitude + st];
endfunction

initServerQ = 500;
initClientQ = 0;
num_clients = 5;
num_servers = 3;
qsz_exponent = 1;
os = 0;

init = vertcat(repmat([initServerQ], num_servers, 1), 
			   repmat([initClientQ], num_clients, 1));
lags = repmat([0.5], 1, num_servers);

hist_mat = vertcat(repmat(initServerQ, num_servers, num_servers),
				   repmat(initClientQ, num_clients, num_servers));

res = ode23d(@q, T, init, lags, hist_mat, num_clients, num_servers, lags, qsz_exponent, os);

% XXX: Add legend
subplot (2, 1, 1);
plot(res.x, res.y(:, 1:num_servers), 'LineWidth', 2);
subplot (2, 1, 2);
plot(res.x,  res.y(:, num_servers + 1 :num_servers + num_clients), 'LineWidth', 2);