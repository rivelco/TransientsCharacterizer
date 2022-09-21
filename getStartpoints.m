function startsVec = getStartpoints(neuron)
    startsVec = {};
    for i = 1:length(neuron)
        prev = i-1;
        if i == 1
            prev = 1;
        end
        if neuron(i) == 1 && neuron(prev) == 0
            startsVec = [startsVec, i];
        end
    end
    startsVec = cell2mat(startsVec);
end
