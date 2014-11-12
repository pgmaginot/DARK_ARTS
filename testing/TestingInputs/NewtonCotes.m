function [x ,w] = NewtonCotes(n)
% n points
% space from -1 to 1
% Newton-Cotes formulas taken from Abramowitz book (publish by US Dept of
% Commerce)
% w must sum to 2
x = linspace(-1,1,n)';
% switch to get fraction of w out for paper
% w = sym( 'w')
% useful w
w = zeros(n,1);
switch n
    case 2
        w(1,1) = 1;
        w(2,1) = 1;
    case 3       
        w(1,1) = 1;
        w(2,1) = 4;  
        w(3,1) = 1;
    case 4        
        w(1,1) = 3;
        w(2,1) = 9;  
        w(3,1) = 9;
        w(4,1) = 3;
    case 5       
        w(1,1) = 7;
        w(2,1) = 32;  
        w(3,1) = 12;
        w(4,1) = 32;  
        w(5,1) = 7; 
    case 6        
        w(1,1) = 19;
        w(2,1) = 75;  
        w(3,1) = 50;
        w(4,1) = 50;  
        w(5,1) = 75;
        w(6,1) = 19;  
    case 7      
        w(1,1) = 41;
        w(2,1) = 216;  
        w(3,1) = 27;
        w(4,1) = 272;  
        w(5,1) = 27;
        w(6,1) = 216;  
        w(7,1) = 41;
    case 8
        w(1,1) = 751;
        w(2,1) = 3577;  
        w(3,1) = 1323;
        w(4,1) = 2989;  
        w(5,1) = 2989;
        w(6,1) = 1323;  
        w(7,1) = 3577;
        w(8,1) = 751;
    case 9
        w(1,1) = 989;
        w(2,1) = 5888;  
        w(3,1) = -928;
        w(4,1) = 10496;  
        w(5,1) = -4540;
        w(6,1) = 10496;  
        w(7,1) = -928;
        w(8,1) = 5888;
        w(9,1) = 989;        
    case 10
        w(1,1) = 2857;
        w(2,1) = 15741;  
        w(3,1) = 1080;
        w(4,1) = 19344;  
        w(5,1) = 5778;
        w(6,1) = 5778;  
        w(7,1) = 19344;
        w(8,1) = 1080;
        w(9,1) = 15741;
        w(10,1) = 2857;  
    case 11
        w(1,1) = 16067;
        w(2,1) = 106300;  
        w(3,1) = -48525;
        w(4,1) = 272400;  
        w(5,1) = -260550;
        w(6,1) = 427368;  
        w(7,1) = -260550;
        w(8,1) = 272400;
        w(9,1) = -48525;
        w(10,1) = 106300;
        w(11,1) = 16067;
    otherwise
        error('Newton-Cotes formula not coded')
end

% normalize weigths

w = 2/sum(w).*w;

return
end