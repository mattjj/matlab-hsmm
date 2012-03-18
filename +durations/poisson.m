classdef poisson < handle
   properties
       % hyperparameters
       k
       theta
       % mean parameter
       lmbda
   end
   
   methods
       function self = poisson(k,theta,lmbda)
           self.k = k;
           self.theta = theta;
           if ~exist('lmbda','var')
               self.resample([]);
           else
               self.lmbda = lmbda;
           end
       end
       
       function resample(self,data)
           self.lmbda = gamrnd(self.k + sum(sum(data-1)),self.theta / (self.theta * size(data,2) + 1));
       end
       
       function val = log_pmf(self,x)
           val = -self.lmbda + (x-1.)*log(self.lmbda) - gammaln(x);
           val(x<1) = log(0.);
       end
       
       function val = pmf(self,x)
           val = poisspdf(x-1,self.lmbda);
       end
       
       function val = log_sf(self,x)
           % this could be coded better, but there's no log incomplete
           % gamma function routine
           val = log(1-poisscdf(x-1,self.lmbda)); 
       end
       
       function sample = rvs(self,num)
           if nargin < 2, num = 1; end
           sample = poissrnd(self.lmbda,1,num);
       end
   end
    
end