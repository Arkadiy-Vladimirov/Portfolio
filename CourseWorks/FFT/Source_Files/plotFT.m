function info = plotFT(hFigure, fHandle, fFTHandle, step,...
                                                       inpLimVec, varargin)
%PLOTFT plots fast Fourier transform of given function (and analytical
%transform, if provided)
%   Detailed explanation goes here
%hFigure - handle of figure to plot on
%fHandle - transformed function handle
%fFTHandle - handle of analitically calculated transform or []
%step - sampling increment
%inpLimVec - window [a b]
%varargin = outLimVec (limits of lamda axis of FFT plot) or NaN
%info - structure containing FFT results

%tiny argument validation block
if (nargin == 6)
    outLimVec = varargin{1};
elseif (nargin ~= 5)
    error('Inadequate number of input values');
else
   outLimVec = [-10 10]; %not calculated yet
end
%end of it

    %calculating FFT
    a = inpLimVec(1); b  = inpLimVec(2);
    T = b - a;
    N = fix(T/step);
    dt = T/N; 
    
    a0 = 0; b0 = b-a; 
    tn = a0:dt:(b0-dt/2);
    dl = 2*pi/T; ln = 0:dl:2*pi*(N-1)/T;
    f0 = @(t) fHandle(t-a);
    f0n = f0(tn);
    F0n = 2*pi/N * fft(f0n);
    Fn = exp(1i*ln*a).*F0n;
    
    if (outLimVec ~= 0)
        c = outLimVec(1); d = outLimVec(2);
        cn = round(c/dl); dn = round(d/dl);
        nN = cn:1:dn;
        ext = nN*dl;
        Fext = Fn(mod(nN,N)+1);
    end
   
    SPlotInfo = get(hFigure,'UserData');
    if ~isempty(SPlotInfo)
        nexttile(1); axIm = gca;
        nexttile(2); axRe = gca;
        outLimVec = SPlotInfo.outLimVec;
    else
        clf(hFigure,'reset');
        hFigure.Name = 'Fourier transform';
        tiledlayout(hFigure,2,1);
        
        axIm = nexttile;
        axIm.XAxisLocation = 'origin'; axIm.YAxisLocation = 'origin';
        axIm.XLabel.Interpreter = 'latex';
        axIm.YLabel.Interpreter = 'latex';
        axIm.XLabel.String = '$\lambda$';
        axIm.YLabel.String = '$Im(F(\lambda))$';
        axIm.YLabel.Rotation = 0;
        axIm.XLim = outLimVec;
        hold on;
        
        axRe = nexttile;
        axRe.XAxisLocation = 'origin'; axRe.YAxisLocation = 'origin';
        axRe.XLabel.Interpreter = 'latex'; 
        axRe.YLabel.Interpreter = 'latex';
        axRe.XLabel.String = '$\lambda$';
        axRe.YLabel.String = '$Re(F(\lambda))$';
        axRe.YLabel.Rotation = 0;
        axRe.XLim = outLimVec;
        hold on;
    end
    
    %plotting
    hNImPl = plot(axIm,ext,imag(Fext));
    legend(axIm,'numerical FT');
    hNRePl = plot(axRe,ext,real(Fext));
    legend(axRe,'numerical FT');
    
    info.AnalyticalSolutionProvided = 'No'; %some info to return
    if ~isempty(fFTHandle)
        hAImPl = fplot(axIm,@(l) imag(fFTHandle(l)),outLimVec);
        legend(axIm,[hNImPl hAImPl],'numerical FT','analytical FT');
        hARePl = fplot(axRe,@(l) real(fFTHandle(l)),outLimVec);
        legend(axRe,[hNRePl hARePl],'numerical FT', 'analytical FT');
        SPlotInfo.hAnPlot.Im = hAImPl; %saving metadata
        SPlotInfo.hAnPlot.Re = hARePl; %saving metadata
        info.AnalyticalSolutionProvided = 'Yes'; %some info to return
    end
    
    %saving metadata
    SPlotInfo.hImAx = axIm;
    SPlotInfo.hReAx = axRe;
    SPlotInfo.hNumPlot.Im = hNImPl;
    SPlotInfo.hNumPlot.Re = hNRePl;
    SPlotInfo.outLimVec = outLimVec;
    set(hFigure,'UserData',SPlotInfo); 
    
    %returning info
    info.nPoints = N;
    info.step = dt;
    info.inpLimVec = inpLimVec;
    info.outLimVec = outLimVec;
end

