classdef gaussian < handle
    properties
        % hyperparameters as in Gelman et al.'s notation
        mu_0 
        lmbda_0
        nu_0
        kappa_0
        % parameters are mean and covariance matrix
        mu
        sigma
        % convenient
        ndim
    end
    
    methods
        function self = gaussian(mu_0,lmbda_0,nu_0,kappa_0,mu,sigma)
            self.mu_0 = mu_0;
            self.lmbda_0 = lmbda_0;
            self.nu_0 = nu_0;
            self.kappa_0 = kappa_0;
            self.ndim = size(lmbda_0,1);
            if ~exist('mu','var') || ~exist('sigma','var')
                self.resample([]);
            else
                self.mu = mu;
                self.sigma = sigma;
            end
        end
        
        function resample(self, data)
            % data is ndim x nsamples
            if isempty(data)
                [self.mu, self.sigma] = self.sample_niw(self.mu_0,self.lmbda_0,self.nu_0,self.kappa_0);
            else
                n = size(data,2);
                % calculate sufficient statistics
                xbar = mean(data,2);
                centered = bsxfun(@minus,data,xbar);
                sumsq = centered*centered';
                % form posterior hyperparameters
                mu_n = self.kappa_0 / (self.kappa_0 + n) * self.mu_0 + n / (self.kappa_0 + n) * xbar;
                kappa_n = self.kappa_0 + n;
                nu_n = self.nu_0 + n;
                lmbda_n = self.lmbda_0 + sumsq + self.kappa_0 * n / (self.kappa_0+n) * (xbar - self.mu_0)*(xbar - self.mu_0)';
                [self.mu, self.sigma] = self.sample_niw(mu_n, lmbda_n, nu_n, kappa_n);
            end
        end
        
        function val = log_likelihood(self,x) 
            x = bsxfun(@minus,x,self.mu);
            val = -1./2. * sum(x.*(self.sigma\x)) - log((2*pi)^(self.ndim/2)*sqrt(det(self.sigma)));
        end
        
        function sample = rvs(self,num)
            sample = bsxfun(@plus,self.mu,cholcov(self.sigma)'*randn(self.ndim,num));
        end

        function plot(self,color,ax)
            if ~exist('color','var'), color = 'b'; end
            if ~exist('ax','var'), ax = gca(); end

            if self.ndim == 2
                util.plot_gaussian_2D(self.mu,self.sigma,color,1,ax);
            else
                error('Plotting not implemented for Gaussians in >2 dimensions.')
            end
        end
    end
    
    methods(Static)
        function [mu, lmbda] = sample_niw(mu_0,lmbda_0,nu_0,kappa_0)
            lmbda = iwishrnd(lmbda_0,nu_0);
            mu = cholcov(lmbda/kappa_0)'*(randn(size(lmbda_0,1),1))+mu_0;
        end
    end
end
