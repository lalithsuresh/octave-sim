clear;
pkg load odepkg;

Tstart = 1;
Tend = 50;
NumT = Tend;
% T = linspace(Tstart, Tend, NumT * 5);
T = [Tstart,Tend];

function qdot = q(t, x, xd,  num_clients, num_servers, lags, qsz_exponent, os, sending_rate)
	qdot = zeros(num_servers + num_clients, 1);
	arrival_rate = ArrivalRate(t, num_clients);
	flow_matrix = zeros(num_clients, num_servers);

	if(t > 11)
		sending_rate(1) = 56/5;
	endif
	for client = 1:num_clients

		% 
		% Flows are divided in proportion to the inverse of their
		% queue-size x service-time products.
		% 
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

		% 
		% Demand allocation begins.
		% 
		total_to_allocate = arrival_rate(client);
		is_full = zeros(num_servers, 1);
		
		% When in backpressure mode, fire away at full possible
		% sending rate.
		% XXX: this check has potential for numerical errors.
		% ... fix later.
		if (x(num_servers + client) > 0)
			total_to_allocate = sum(sending_rate);
		endif

		incoming_demand = total_to_allocate;
		% The inner loop allocates the current demand (which is either the
		% arrival rate when not in backpressure, or the full sending rate when
		% in backpressure) across all servers. At the end of each iteration,
		% any allocation which exceeds the sending rate for a server is
		% accumulated as the backlogRate. The outer loop repeats this
		% iteration either until all sending rates have been saturated or
		% until the demand has been fully allocated across servers.
		for k = 1:num_servers
			backlogRate = 0;
			WeightToUse = TotalWeight;

			for server = 1:num_servers
				if (is_full(server) == 0)
					flow_fraction = qmu_inverse(server)/(WeightToUse);
					flow_matrix(client, server) = flow_matrix(client, server) + total_to_allocate * flow_fraction;

					if(flow_matrix(client, server) > sending_rate(server))
						backlogRate = backlogRate + flow_matrix(client, server) - sending_rate(server);
						flow_matrix(client, server) = sending_rate(server);
						is_full(server) = 1;
						TotalWeight = TotalWeight - qmu_inverse(server);
					endif
				endif
			endfor

			% If there is a residual demand, that is allocated in the next
			% iteration.
			if (backlogRate != 0)
				total_to_allocate = backlogRate;
			else
				break
			endif
		endfor

		if(backlogRate > 0)
			t, backlogRate, flow_matrix
		endif

		% Update backlog size at client
		total_drain = sum(flow_matrix(client, :));
		% backlogRate, total_drain, incoming_demand 
		% assert(backlogRate + total_drain <= incoming_demand + 0.1 && backlogRate + total_drain >= incoming_demand - 0.1);

		% if (x(num_servers + client) <= 0)
		% 	total_drain = 0;
		% endif
		qdot(num_servers + client) = arrival_rate(client) - total_drain;
	endfor

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

	if(t > 11)
		ArrivalRate = repmat(lamb, 1, num_clients);
	elseif(t > 10)
		ArrivalRate = repmat(30, 1, num_clients);
	endif

endfunction

% xxx: use sinusoid
function ServiceRate = ServiceRate(t, num_servers)
	
	ServiceRate = repmat(50, 1, num_servers);

	if(t > 11)
		ServiceRate = [56;50;50];
		% ServiceRate = [30;70;50];
	elseif(t > 10)
		% ServiceRate = [30;70;50];
		ServiceRate = [50;50;50];
	endif

	% st = 50;
	% amplitude = 5;
	% ServiceRate = [sin(t/10 + pi/4) * amplitude + st; 
	% 			   sin(t/10 + pi/2) * amplitude + st; 
	% 			   sin(t/10 + 3* pi/4) * amplitude + st];
endfunction

initServerQ = 10;
initClientQ = 0;
num_clients = 5;
num_servers = 3;
qsz_exponent = 3;
os = 0;

init = vertcat(repmat([initServerQ], num_servers, 1), 
			   repmat([initClientQ], num_clients, 1));
lags = repmat([0.5], 1, num_servers);
hist_mat = vertcat(repmat(initServerQ, num_servers, num_servers),
				   repmat(initClientQ, num_clients, num_servers));
sending_rates = repmat([10], 1, num_servers);

options = odeset ('InitialStep', 0.1, 'MaxStep', 0.1, 'RelTol', 1.0e-6, 'AbsTol', 1.0e-6); 
res = ode23d(@q, T, init, lags, hist_mat, options,
			 num_clients, num_servers, lags, qsz_exponent, os, sending_rates);

% XXX: Add legend
subplot (2, 1, 1);
plot(res.x, res.y(:, 1:num_servers), 'LineWidth', 2);
axis ([Tstart Tend 0 150]);
subplot (2, 1, 2);
plot(res.x,  res.y(:, num_servers + 1 :num_servers + num_clients), 'LineWidth', 2);
axis ([Tstart Tend 0 5]);