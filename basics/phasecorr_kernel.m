function ckernel = phasecorr_kernel(kernel)

dim = size(kernel,1);

for k = 1:dim
    ckernel(k,:) = phasecorr(kernel(k,:)); 
end

end