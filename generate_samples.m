function varsamples = generate_samples(Nvar, percent_var, Nsamples)

Nvarsample = round(Nvar.*percent_var./100);
varsamples = zeros(Nvarsample, Nsamples);
for k=1:Nsamples
    varsamples(:,k) = randsample(Nvar, Nvarsample);
end
