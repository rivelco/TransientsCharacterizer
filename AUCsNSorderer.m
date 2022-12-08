function [scores, idxs] = AUCsNSorderer(ensemble, property)
    result = property{ensemble};
    result = result';
    result(isnan(result)) = 0;
    final = sortrows(result, [2, 1], 'descend');
    idxs = final(:, 1);
    scores = final(:, 2);
end