clear, clc

n = 10000;
Actives = zeros(n, 1);
Nona = zeros(n, 1);

tic
for i = 1:n
    [a, b] = VectorsAnalyzer;
    Actives(i) = a;
    Nona(i) = b;
end
toc

X = [Actives; Nona];
grp = [ones(size(Actives)); 2.*ones(size(Nona))];

boxplot(X, grp)
