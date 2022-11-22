function [some] = plotModelData(data, modelUp, modelDown, xUp, xDown)
    figure(1), clf
    a = modelUp(xUp);
    a2 = a/max(a);
    b = modelDown(xDown);
    b2 = b/max(b);
    c = data;
    c2 = c/max(c);
    plot(xUp, a2, 'LineWidth', 2, 'Color', 'red')
    hold on
    plot(xDown, b2, 'LineWidth', 2, 'Color', 'blue')
    hold on
    plot(1:max(xDown), c2, '.', 'Color', 'black', 'MarkerFaceColor','black')
    xline(100, 'LineWidth', 2, 'Color', 'green')
    some = 0;
end