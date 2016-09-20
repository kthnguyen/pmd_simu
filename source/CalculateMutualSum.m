function [output] = CalculateMutualSum(A,B)
%Calculate sum of each element in an 1-row vector A to an 1-row vector B
%Input: 1-row vectors A and B
%Output: Matrix size(A)*size(B)

szA = length(A);
szB = length(B);

repA = repmat(A',1,szB);
repB = repmat(B,szA,1);

output = repA + repB;

end
