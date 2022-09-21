function endsVec = getEndpoints(neuron)
    endsVec = {};
    for i = 1:length(neuron)
        prev = i-1;
        if i == 1
            prev = 1;
        end
        if neuron(i) == 0 && neuron(prev) == 1
            endsVec = [endsVec, i];
        end
    end
    endsVec = cell2mat(endsVec);
end
