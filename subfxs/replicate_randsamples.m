function replicate_randsamples(outpath)

rng('shuffle')
goodsubj = [1:128];    
Nsub   = length(goodsubj);
p = randperm(Nsub);
isuj1 = p(1:Nsub/2);
isuj2 = p(Nsub/2+1:Nsub);
goodsubj1 = sort(goodsubj(isuj1));
goodsubj2 = sort(goodsubj(isuj2));
cd(outpath)
save replicates_goodsubj goodsubj1 goodsubj2
