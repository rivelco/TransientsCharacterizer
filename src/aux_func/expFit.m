function [fitUp, gofUp, fitDown, gofDown] = expFit(data, spkUp, spkDown)
    yUp = data(spkUp);
    xUp = spkUp;
    yDown = data(spkDown);
    xDown = spkDown;

    ft = fittype('exp2');
    optsUp = fitoptions( 'Method', 'NonlinearLeastSquares' );
    optsUp.Display = 'Off';
    optsUp.Lower = [-Inf -Inf -Inf 0];
    optsUp.Upper = [Inf Inf Inf 0];

    optsDown = fitoptions( 'Method', 'NonlinearLeastSquares' );
    optsDown.Display = 'Off';
    optsDown.Lower = [-Inf -Inf -Inf 0];
    optsDown.Upper = [Inf Inf Inf 0];
    
    [xData, yData] = prepareCurveData( xUp, yUp );
    [fitUp, gofUp] = fit(xData, yData, ft, optsUp);
    [xData, yData] = prepareCurveData( xDown, yDown );
    [fitDown, gofDown] = fit(xData, yData, ft, optsDown);
end