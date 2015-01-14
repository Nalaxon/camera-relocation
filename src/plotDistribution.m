function plotDistribution(x, pd)
hold on;
for i=1:1%size(pd,1)
    y = 1/sqrt(2*pi*pd(i,2))*exp(-1/2*((x'-pd(i,1)) / sqrt(pd(i,2)) ).^2);
    mu = pd(i,1)
    sig = pd(i,2)
    plot(x', y);
end
hold off;
end


