function params = best_parameters(X)
	% X: dataset, needs to be an output of Cox_data_processed.m

	c_index_tau = cell(1,100); 

	count = 0;
	count_tau = 0;
	tot = 250000; % number of total experiments

	iterations = cell(1,100);

	% Loop for parameter selection
	range_tau = 0.01:0.01:1;
	for tau = range_tau % Cycle on tau

	    c_index_rho = zeros(50);
	    iterations_rho = zeros(50);
	    count_rho = 0;

	    range_rho = linspace(0,10,50);
	    for rho = range_rho  % Cycle on rho

	        c_index_lambda = zeros(1,50);
	        iterations_lambda = zeros(1,50);
	        count_lambda = 0;

	        range_lambda = linspace(0,10,50);
	        for lambda = range_lambda  % Cycle on lambda
	            
	            % Perform S2GC
	            [wm, funcValMv,G,iter] = S2GC(X, lambda, rho, tau, ViewIdx);
	            
	            % Compute concordance index
	            c = concordance_index(X.sorty , X.X*wm ,~X.cens,1);

	            count_lambda = count_lambda + 1; % Increase lambda counter

	            c_index_lambda(count_lambda) = c;

	            iterations_lambda(count_lambda) = iter;

	            count = count + 1;  % Increase counter
	            fprintf('###### Iteration ######## %d su %d %s \n\n ', ...
                    count, tot, datetime(now,'ConvertFrom','datenum'));
	            
	        end

	        count_rho = count_rho + 1;  % Increase rho counter

	        c_index_rho(count_rho,:) = c_index_lambda;

	        iterations_rho(count_rho,:) = iterations_lambda;

	    end

	    count_tau = count_tau + 1;  % Increase tau counter

	    c_index_tau(count_tau) = {c_index_rho};

	    iterations(count_tau) = {iterations_rho};

	end

	[optimal_c_index_tau, i_best_params_tau] = max(c_index_tau);
	[optimal_c_index_rho, i_best_params_rho] = max(c_index_rho);
	[optimal_c_index_lambda, i_best_params_lambda] = max(c_index_lambda);

	optimal_tau = range_tau(i_best_params_tau);
	optimal_rho = range_rho(i_best_params_rho);
	optimal_lambda = range_lambda(i_best_params_lambda);

	params = struct('tau', optimal_tau, 'rho', optimal_rho, 'lambda', ...
        optimal_lambda);

end