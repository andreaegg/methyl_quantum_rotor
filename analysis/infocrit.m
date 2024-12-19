function metrics = infocrit(V,Vfit,nrparam)

N = length(V);
K = nrparam + 2;
res = norm(V-Vfit);

metrics.aic  = N*log(res.^2/N) + 2*K;
metrics.aicc = N*log(res.^2/N) + 2*K*N/(N-K-1);
metrics.bic  = N*log(res.^2/N) +K*log(N);

end