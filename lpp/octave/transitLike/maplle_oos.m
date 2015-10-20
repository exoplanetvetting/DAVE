function oosX = maplle_oos(X,mapping,dim)
%Code to create an out of sample set given 
%an already known set of mapping with LLE
%mapping comes after running LLE or compute_mapping(X','LLE')

s=size(X);
oosX=zeros(s(1),dim);

for i=1:s(1)
    
    oosX(i,:)=out_of_sample(X(i,:),mapping);
    
end


