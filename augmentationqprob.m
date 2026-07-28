function[prob]=augmentationqprob(m,a,b,coeff,k,delta,delta2,delta3,delta4)
% Runs m times augmentationq for a singular quadratic matrix polynomial
% and provides the empirical probability of getting the right
% number of true evals
% a is the number of true finite evals
% b is the number of true infinite evals

rng(0);
number=0;
for i=1:m
    [~,evals,inf]=augmentationq(coeff,k,delta,delta2,delta3,delta4);
    t1=length(evals);
    t2=length(inf);
    if t1==a && t2==b
        number=number+1;
    end
end
prob=number/m;