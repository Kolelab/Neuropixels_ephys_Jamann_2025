function pwr = bootstrap_power(x,n,test,alpha)
%bootstrap_power. Returns the power when from a dataset n samples are chosen
%
%   pwr = bootstrap_power(x,n,test,alpha)
%
% 2025, Alexander Heimel

if nargin<4 || isempty(alpha)
    alpha = 0.05;
end
if nargin<5 || isempty(test)
    test = 'ttest';
end

switch test
    case 'ttest'
        test = @my_ttest;
    case 'signrank'
        test = @signrank;
end

n_tries = 10000;

pwr = 0;
for i = 1:n_tries
    ind = randi(length(x),[1 n]);
    p = test(x(ind));
    if p<alpha
        pwr = pwr+1;
    end

end % i

pwr = pwr / n_tries;
end

function p = my_ttest(x)
[~,p] = ttest(x);
end