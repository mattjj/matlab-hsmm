function plot_gaussian_2D(mu, lmbda, color, centermarker, ax)
    if ~exist('color','var'), color = 'b'; end
    if ~exist('centermarker','var'), centermarker = 0; end
    if ~exist('ax','var'), ax = gca(); end

    t = [0:0.01:2*pi,0];
    circle = [sin(t);cos(t)];
    ellipse = chol(lmbda)*circle;
    
    if strcmp(class(color),class('string'))
        if centermarker
            plot(ax,mu(1),mu(2),strcat(color,'+','MarkerSize',10));
        end
        plot(ax,ellipse(1,:)+mu(1),ellipse(2,:)+mu(2),strcat(color,'-'),'LineWidth',2);
    elseif strcmp(class(color),class([1,2,3]))
        if centermarker
            plot(ax,mu(1),mu(2),'+','Color',color,'MarkerSize',10);
        end
        plot(ax,ellipse(1,:)+mu(1),ellipse(2,:)+mu(2),'-','Color',color,'LineWidth',2);
    else
        error('Unrecognized color variable.')
    end

end
