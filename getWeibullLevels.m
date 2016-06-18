params = [];
levels = [0.67, 0.75, 0.9];


for n = 1:5
    
    newStimLevel = zeros(1,3);
    for ii = 1:3
        newStimLevel(ii) = results.weibullFitParams(n).x(find(results.weibullFitParams(n).y > levels(ii), 1))
    end
    
    params = [params; newStimLevel]
end
