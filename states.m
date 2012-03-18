classdef states < handle
    properties
        T
        state_dim
        obs_distns
        dur_distns
        transition_distn
        initial_distn
        trunc
        
        stateseq
        durations
        stateseq_norep
    end
    
    methods
        function self = states(T,state_dim,obs_distns,dur_distns,transition_distn,initial_distn,trunc)
            if ~exist('trunc','var'), self.trunc = T+1; else self.trunc = trunc; end
            
            self.T = T;
            self.state_dim = state_dim;
            self.obs_distns = obs_distns;
            self.dur_distns = dur_distns;
            self.transition_distn = transition_distn;
            self.initial_distn = initial_distn;
            
            self.generate();
        end
        
        function resample(self,data)
            % generate duration pmf, sf, and likelihood values
            possible_durations = 1:self.T;
            aDl = zeros(self.T, self.state_dim);
            aDsl = zeros(self.T,self.state_dim);
            aBl = zeros(self.T,self.state_dim);
            for idx=1:self.state_dim
                aDl(:,idx) = self.dur_distns{idx}.log_pmf(possible_durations)';
                aDsl(:,idx) = self.dur_distns{idx}.log_sf(possible_durations)';
                aBl(:,idx) = self.obs_distns{idx}.log_likelihood(data)';
            end
            % run backwards message passing
            [betal, betastarl] = self.messages_backwards(log(self.transition_distn.A),aBl',aDl',aDsl',self.trunc);
            % sample forwards
            self.sample_forwards(aBl',betal,betastarl);
        end
        
        function sample_forwards(self,aBl,betal,betastarl)
            T = size(aBl,2);
            
            stateseq = zeros(1,T);
            durations = [];
            stateseq_norep = [];
            
            idx = 1;
            A = self.transition_distn.A;
            nextstate_unsmoothed = self.initial_distn.pi_0;
            
            while idx <= self.T
                logdomain = betastarl(:,idx) - max(betastarl(:,idx));
                nextstate_distr = exp(logdomain) .* nextstate_unsmoothed;
                % due to numerical issues, might be all zero.
                % if that's the case, there's no easy fix, so we'll just
                % follow the messages.
                if all(nextstate_distr == 0)
                    nextstate_distr = exp(logdomain);
                end
                state = util.sample_discrete(nextstate_distr);
                
                durprob = rand();
                dur = 0; % always incremented at least once
                breaking = 0; % see section at line 88
                while durprob > 0 && ~breaking
                    % in the body of this loop, we're really considering
                    % dur+1 being the duration for this state.
                    p_d_marg = self.dur_distns{state}.pmf(dur+1);
                    if p_d_marg == 0
                        dur = dur + 1;
                        continue;
                    end
                    if idx+dur <= self.T
                        p_d = exp(sum(aBl(state,idx:idx+dur),2) + betal(state,idx+dur) - betastarl(state,idx)) * p_d_marg;
                        durprob = durprob - p_d;
                        dur = dur + 1;
                    else
                        % we're out of data, so we need to sample a
                        % duration conditioned on us having lasted at least
                        % this long. the likelihood contributes the same to
                        % all possibilities, so we can just sample from the
                        % prior (conditioned on it being at least this
                        % long). rejection sampling is easiest!
                        while 1
                            lastdur = self.dur_distns{state}.rvs();
                            if lastdur > dur
                                dur = lastdur;
                                breaking = 1;
                                break;
                            end
                        end
                    end
                end
                
                stateseq(idx:min(idx+dur-1,T)) = state;
                stateseq_norep = [stateseq_norep, state];
                durations = [durations, dur];
                
                nextstate_unsmoothed = A(state,:)';
                
                idx = idx + dur;
            end
            
            
            self.stateseq = stateseq;
            self.durations = durations;
            self.stateseq_norep = stateseq_norep;
        end
        
        function obs = generate(self)
            % generate states
            idx = 1;
            nextstate_distr = self.initial_distn.pi_0;
            A = self.transition_distn.A;
            
            stateseq = [];
            stateseq_norep = [];
            durations = [];
            
            while idx <= self.T
                % sample a state
                state = util.sample_discrete(nextstate_distr);
                % sample a druation for that state
                duration = self.dur_distns{state}.rvs();
                % save everything
                stateseq_norep = [stateseq_norep, state];
                durations = [durations, duration];
                stateseq = [stateseq, state * ones(1,duration)];
                % set up next state distribution
                nextstate_distr = A(state,:)';
                % update index
                idx = idx + duration;
            end
            
            self.stateseq_norep = stateseq_norep;
            self.durations = durations;
            self.stateseq = stateseq(1:self.T);
            
            % generate observations
            % not pre-allocated just because we don't know the obs dimension here
            obs = [];
            for t=1:self.T
                obs = [obs, self.obs_distns{self.stateseq(t)}.rvs(1)];
            end
        end
    end

    methods(Static)
        function [betal, betastarl] = messages_backwards(Al,aBl,aDl,aDsl,trunc)
            T = size(aBl,2);
            state_dim = size(Al,1);

            betal = zeros(state_dim,T);
            betastarl = zeros(state_dim,T); 
            
            for t=T:-1:2
                betastarl(:,t) = util.logsumexp(betal(:,t:min(T,t+trunc-1)) + cumsum(aBl(:,t:min(T,t+trunc-1)),2) + aDl(:,1:min(trunc-1,T-t+1)));
                if T-t+1 < trunc
                    betastarl(:,t) = util.logsumexp([betastarl(:,t),sum(aBl(:,t:end),2) + aDsl(:,T-t+1)]); % censoring calc
                end
                betal(:,t-1) = util.logsumexp(bsxfun(@plus,betastarl(:,t)',Al));
            end
            t = 1;
            betastarl(:,t) = util.logsumexp(betal(:,t:min(T,t+trunc-1)) + cumsum(aBl(:,t:min(T,t+trunc-1)),2) + aDl(:,1:min(trunc-1,T-t+1)));
            if T-t+1 < trunc
                betastarl(:,t) = util.logsumexp([betastarl(:,t),sum(aBl(:,t:end),2) + aDsl(:,T-t+1)]); % censoring calc
            end
        end
    end
end
