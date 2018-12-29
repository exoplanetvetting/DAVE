
function sum = testFunc(mapStructFile)
%A test function that does some kind of computation

fprintf(1, mapStructFile)
mapStruct = load(mapStructFile);

sum = 0
for i = 1:1e1
    for j = 1:1e2
        sum = sum + 1
    end

end    

end