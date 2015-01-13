function plotDistribution(x, pd, H)
hold on;
for i=1:1%size(pd,1)
    y = 1/sqrt(2*pi*pd(2*i)^2)*exp(-(x'-pd(2*i-1)).^2/(2*pd(2*i-1)));
    plot(x', y);
end
hold off;
end


