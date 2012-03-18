classdef transitions < handle
    properties
        % Dirichlet hyperparameters
        alpha
        gamma
        state_dim
        % beta (mean transition row) and transition matrix parameters
        beta
        fullA
        A
    end
    
    methods
        function self = transitions(state_dim, alpha, gamma, beta, A, fullA)
            self.state_dim = state_dim;
            self.alpha = alpha;
            self.gamma = gamma;
            if ~exist('beta','var') || ~exist('A','var') || ~exist('fullA','var')
                self.resample([]);
            else
                self.beta = beta;
                self.A = A;
                self.fullA = fullA;
            end
        end
        
        function resample(self,states_norep)
            if size(states_norep,2) < 2
                % there is no transition data, just sample forward
                self.beta = gamrnd(ones(1,self.state_dim)*self.alpha/self.state_dim,1);
                self.beta = self.beta / sum(sum(self.beta));
                
                self.fullA = gamrnd(self.gamma * self.beta(ones(self.state_dim,1),:),1);
                self.A = (1-eye(self.state_dim)) .* self.fullA;
                self.fullA = bsxfun(@rdivide,self.fullA,sum(self.fullA,2));
                self.A = bsxfun(@rdivide,self.A,sum(self.A,2));
            else
                % make 2D array of transition counts
                data = zeros(self.state_dim,self.state_dim);
                for idx=1:size(states_norep,2)-1
                    data(states_norep(idx),states_norep(idx+1)) = data(states_norep(idx),states_norep(idx+1)) + 1;
                end
                % complete the data using instantiated fullA parameters
                froms = sum(data,2);
                self_transitions = zeros(1,length(froms));
                for idx=1:length(froms)
                    self_transitions(idx) = sum(geornd(1-self.fullA(idx,idx),froms(idx),1));
                end
                augmented_data = data + diag(self_transitions);
                % then, compute the m auxiliary variables (as described in 
                % the E.B. Fox thesis)
                m = zeros(self.state_dim,self.state_dim);
                for rowidx=1:self.state_dim
                    for colidx=1:self.state_dim
                        for idx=1:augmented_data(rowidx,colidx)
                            m(rowidx,colidx) = m(rowidx,colidx) + (rand() < self.alpha * self.beta(colidx) / ((idx-1)+self.alpha*self.beta(colidx)));
                        end
                    end
                end
                % resample mother (beta)
                self.beta = gamrnd(self.alpha / self.state_dim + sum(m),1);
                self.beta = self.beta / sum(sum(self.beta));
                % resmaple children (fullA and A)
                self.fullA = gamrnd(bsxfun(@plus,self.gamma * self.beta,augmented_data),1);
                self.fullA = bsxfun(@rdivide,self.fullA,sum(self.fullA,2));
                self.A = self.fullA .* (1-eye(self.state_dim));
                self.A = bsxfun(@rdivide,self.A,sum(self.A,2));
            end
        end
    end
end
