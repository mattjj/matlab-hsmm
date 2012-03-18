%% set up parameters
% prior stuff for synthetic data
N = 4;
T = 500;
obs_dim = 2;
durparams = [12,18,25,30];
dur_hypparams = {8,5};
obs_hypparams = {zeros(obs_dim,1),eye(obs_dim),obs_dim+3,0.1};

% inference
Nmax = 10; % weak limit approximation parameter
Niter = 100;
plot_results = 1;
plot_results_every = 10;

%% create truth model
for state=1:N
    truth_obs_distns{state} = observations.gaussian(obs_hypparams{:});
    truth_dur_distns{state} = durations.poisson(dur_hypparams{:},durparams(state));
end

truthmodel = hsmm(T,truth_obs_distns,truth_dur_distns);

%% generate data from truth model
[data, labels] = truthmodel.generate();

if plot_results
    fig = figure();
    truthmodel.plot(data,fig);
    title('TRUTH');
end

%% create posterior model
for state=1:Nmax
    obs_distns{state} = observations.gaussian(obs_hypparams{:});
    dur_distns{state} = durations.poisson(dur_hypparams{:});
end

posteriormodel = hsmm(T,obs_distns,dur_distns);

%% do posterior inference
if plot_results
    fig = figure();
end

for iter=1:Niter
    posteriormodel.resample(data);
    util.print_dot(iter,Niter);

   if plot_results && (mod(iter,plot_results_every) == 0)
       posteriormodel.plot(data,fig);
       title(sprintf('SAMPLED at iteration %d',iter));
       drawnow
   end
end


