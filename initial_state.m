classdef initial_state < handle
    properties
        state_dim
        rho
        pi_0
    end
    
    methods
        function self = initial_state(state_dim,rho,pi_0)
            self.state_dim = state_dim;
            self.rho = rho;
            
            if ~exist('pi_0','var')
                self.resample([]);
            else
                self.pi_0 = pi_0;
            end
        end
        
        function resample(self,initial_states)
            data = full(sum(sparse(initial_states,1:length(initial_states),1,...
                self.state_dim,length(initial_states),length(initial_states)),2));
            self.pi_0 = gamrnd(self.rho / self.state_dim + data,1);
            self.pi_0 = self.pi_0 / sum(self.pi_0);
        end
    end
end