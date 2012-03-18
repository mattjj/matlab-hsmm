classdef hsmm < handle
% convenient wrapper for all the parts of the hsmm definition, so there's only
% one resample() method to call

    properties
        T
        state_dim
        obs_distns
        dur_distns
        trans_distn
        init_state_distn
        states
    end

    methods
        function self = hsmm(T,obs_distns,dur_distns,trunc,trans_distn,init_state_distn)
            state_dim = length(obs_distns);
            self.state_dim = state_dim;
            self.T = T;

            self.obs_distns = obs_distns;
            self.dur_distns = dur_distns;
            
            if ~exist('trans_distn','var')
                self.trans_distn = transitions(state_dim,4,8);
            else
                self.trans_distn = trans_distn;
            end

            if ~exist('init_state_distn','var')
                self.init_state_distn = initial_state(state_dim,2);
            else
                self.init_state_distn = init_state_distn;
            end

            if ~exist('trunc','var')
                self.states = states(T,state_dim,obs_distns,dur_distns,self.trans_distn,self.init_state_distn);
            else
                self.states = states(T,state_dim,obs_distns,dur_distns,self.trans_distn,self.init_state_distn,trunc);
            end
        end

        function resample(self,obs)
            % resample obsparams
            for state=1:length(self.obs_distns)
                self.obs_distns{state}.resample(obs(:,self.states.stateseq == state));
            end

            % resample durparams
            for state=1:length(self.dur_distns)
                self.dur_distns{state}.resample(self.states.durations(self.states.stateseq_norep == state));
            end

            % resample transitions
            self.trans_distn.resample(self.states.stateseq_norep);

            % resample pi_0
            self.init_state_distn.resample(self.states.stateseq(1));

            % resample states
            self.states.resample(obs);
        end

        function [obs, stateseq] = generate(self)
            obs = self.states.generate();
            stateseq = self.states.stateseq;
        end

        function plot(self,data,fig)
            % plots 2D projection of data (if given) and obs distns (if
            % implemented)           
            if ~exist('fig','var'), fig = gcf(); end

            figure(fig);
            clf();
            
            colors = jet(self.state_dim);
            colormap(colors);

            % plot data and obs_distns
            subplot(3,1,[1,2]);
            ax = gca();
            set(ax,'NextPlot','add');
            
            if exist('data','var')
                plot(ax,data(1,:),data(2,:),'k.'); hold on;
            end
            used = zeros(1,self.T);
            used(self.states.stateseq) = 1;
            for state=1:self.state_dim
                if used(state)
                    self.obs_distns{state}.plot(colors(state,:),ax);
                end
            end
            
            % plot state sequence
            subplot(3,1,3);
            image(self.states.stateseq);
            xlabel('sequence index');
            
            % set top axis as current axis, so that things like title can
            % be set easily by caller
            subplot(3,1,[1,2]);
        end
    end

end
