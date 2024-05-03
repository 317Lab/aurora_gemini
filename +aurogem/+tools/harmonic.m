% coefficients provided by A. Mule, October 2023
function f = harmonic(X,Y,pars,order)
arguments
    X (:,:) double {mustBeNumeric}
    Y (:,:) double {mustBeNumeric}
    pars (1,:) double {mustBeNumeric}
    order (1,1) int32 {mustBePositive}
end

num_pars = numel(pars);
assert(num_pars >= 2*order,'pars lenght must equal twice the order.')
assert(mod(num_pars,2)==0,'pars must have even length.')
assert(isequal(size(X),size(Y)),'X2 and X3 must have equal sizes.')

f = zeros(size(X));
order = double(order);
for n = 1:order
    for i = 0:n
        if mod(i,2)==0
            c_even = factorial(n)*(-1)^(i/2) / (factorial(i)*factorial(n-i));
            f = f + pars(2*n-1)*c_even*X.^i.*Y.^(n-i);
        else
            c_odd = factorial(n)*(-1)^((i-1)/2) / (factorial(i)*factorial(n-i));
            f = f + pars(2*n)*c_odd*X.^i.*Y.^(n-i);
        end
    end
end

end