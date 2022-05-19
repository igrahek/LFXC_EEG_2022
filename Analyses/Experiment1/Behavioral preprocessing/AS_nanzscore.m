function [newmat] = AS_nanzscore(mat)


% if n-by-1 array, transposes and then zscores:
if size(mat,1)==1
    mat = mat';
    needsTranspose = 1;
else
    needsTranspose = 0;
end
  
for ci =1:size(mat,2)
    newmat(:,ci) = (mat(:,ci)-nanmean(mat(:,ci)))./nanstd(mat(:,ci));
end

% if n-by-1 array, transposes back to original dimensions:
if needsTranspose
    mat = mat';
end