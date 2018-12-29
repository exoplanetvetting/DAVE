function [ ] = print_tmetricResults( decision,d,label)
%Print the results of the LPP test
%Run lpp_similar_testsample or something that gives a decison for each
%object.
%d is an integer type for the objects.

opts=[0,1,2,3,4,5,6,7,8,9,10];

disp(label)
disp('Train Set   -----   Test Set ')
disp('--Percent Not Transit Like--')

for i=1:length(opts)

    per=100*length(decision(decision' & d==opts(i)))/length(decision(d==opts(i)));
    n=length(decision(d==opts(i)));
    
    %And the same for the test set  which have opts+100
    to=opts(i)+100;
    testper=100*length(decision(decision' & d==to))/length(decision(d==to));
    testn=length(decision(d==to));
    
    s=sprintf('%i | %4.1f   | %i  ||  %4.1f  | %i \n',opts(i),per,n,testper,testn);
    disp(s)

end
disp('----------')



end

