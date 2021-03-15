function metric = rmsd(A,B)
metric = sqrt((1/numel(A))*norm(A-B));
end