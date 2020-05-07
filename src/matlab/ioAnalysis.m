%% Analyse output from NEURON
%   This file is very un-MATLAB-esque.
%       : there is one main function (named after the file),
%       : but lots of other functions within the file too.
%       : MATLAB prefers 1 function per file.
%       : I don't mind a monolithic approach, but unfortunately can't split
%       : functions to analyse and the data run functions due to MATLAB's
%       : scope restrictions.
%       : The easiest way to deal with scrolling is 'ctrl'/'cmd' + '=' to 
%       : collapse everything and 'ctrl'/'cmd' + 'shift' + '=' to expand.

%% Suppress errors
%#ok<*AGROW>    % grow in loop, consider preallocation 
%           --> fast enough
%#ok<*NASGU>    % might be unused 
%           --> not every param needed, but useful when duplicating
%           templates
%#ok<*DEFNU>    % function might be unused
%           --> not every function is run, especially analysis functions
%#ok<*ASGLU>    % return parameter unused 
%           --> best to know potentially returned params for on-demand use
%           --> little cost involved

%% MAIN
function [result] = ioAnalysis()
    global quickPlot toPlot toSave folderName overWrite
    close all;
    %% settings
    % Change default axes fonts.
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultAxesFontSize', 18)

    % Change default text fonts.
    set(0,'defaultTextFontName', 'Arial')
    set(0,'defaultTextFontSize', 18)

    % The properties we've been using in the figures
    set(0,'defaultLineLineWidth',2);   % set the default line width
    set(0,'defaultLineMarkerSize',5); % set the default line marker size
    set(0,'defaultAxesLineWidth',2);   % set the default line width
    set(groot, 'defaultAxesTickDir', 'out'); % The only other option is 'in'
    set(groot,  'defaultAxesTickDirMode', 'manual');
    % Set the default Size for display
    defpos = get(0,'defaultFigurePosition');
    %set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

    % Set the defaults for saving/printing to a file
    set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
    set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
    set(0, 'defaultFigureUnits', 'Inches', ...
        'defaultFigurePosition', [[0,0], [8.27,11.69].*[1,0.5]],...
        'defaultFigurePaperUnits', 'Inches', ...
        'defaultFigurePaperSize', [8.27,11.69])
    
    quickPlot = 0;  %# change to 0 to save as png and eps
    toPlot = 1;
    toSave = 0;
    overWrite = 0;  %# should image files be overwritten? 
                    %# (sometimes off if there is a possibility that 
                    %# different figures have the same saved name)
    folderName='';
    if(toSave)
        toPlot=1; 
    end
    result = 1;

    distal_balanced_synapses_figure();
    proximal_balanced_synapses_figure();
    Data_Increasing_synapse_numbers_proximal();
    Data_Increasing_synapse_numbers_distal();
    result = n_index_figure();
    
end

%% OPTIONS FOR VARIABLES
function [split,name,trimBegin,trimEnd] = getSplitOnNameInPt()
    %split on values between name and inhibition
    split='_[\[()\],e\-\.\d^_]+_InPt';
    name='';
    trimBegin = 1;
    trimEnd = 5;
end
function [split,name,trimBegin,trimEnd] = getSplitOnKCC2()
    %split on KCC2
    split='KCC2';
    name='KCC2';
    trimBegin = 0;
    trimEnd = 0;
end
function [split,name,trimBegin,trimEnd] = getSplitOnNoise()
    %split on values of noise
    split = 'noise([\d\.]*,[\d\.]*';
    name='Noise';
    trimBegin = 6;
    trimEnd = 0;
end
function [split,name,trimBegin,trimEnd] = getSplitOnWeights()
    %split on values of noise
    split = 'weights(-*[\d\.]*,-*[\d\.]*';
    name='Weights';
    trimBegin = 8;
    trimEnd = 0;
end
function [degree,name,trimBegin,trimEnd] = getDegreeOfInhibition()
    %degree of inhibition strength
    
%     degree='InPt\d+_';      
    degree = 'InPt[.\d]+[(_]';
    name='Inhibition Strength';
    trimBegin = 4;
    trimEnd = 1;
end
function [degree,name,trimBegin,trimEnd] = getDegreeOfNameInPt()
    %degree of values between name and inhibition
    degree='_[^_.]+_InPt';
    name='';
    trimBegin = 1;
    trimEnd = 5;
end
function [keep, keepname, funcname] = getKeepProximal()
    %# proximal in general
    keep='proximal*';     
    keepname='proximal';
    funcname='getKeepProximal';
end
function [keep, keepname, funcname] = getKeepproximal()
    %# proximal in general
    keep='proximal*';
    keepname='proximal';
    funcname='getKeepproximal';
end
function [keep, keepname, funcname] = getKeepproximalKCC2()
    keep='proximal_KCC2*';            % proximal with KCC2
    keepname='proximal KCC2';
    funcname='getKeepproximalKCC2';
end
function [keep, keepname, funcname] = getKeepproximalNoKCC2()
    keep='proximal_[^K][^C]+[^2][^_].*';   % proximal not KCC2
    keepname='proximal no KCC2';
    funcname='getKeepproximalNoKCC2';
end
function [keep, keepname, funcname] = getKeepDistal()
    %# distal in general
    keep='distal*';     
    keepname='distal';
    funcname='getKeepDistal';
end
function [keep, keepname, funcname] = getKeepDistalKCC2()
    keep='distal_KCC2*';          
    keepname='distal KCC2';
    funcname='getKeepDistalKCC2';
end
function [keep, keepname, funcname] = getKeepDistalNoKCC2()
    keep='distal_[^K][^C]+[^2][^_].*';
    keepname='distal no KCC2';
    funcname='getKeepDistalNoKCC2';
end
function [keep, keepname, funcname] = getKeepSomatic()
    %# proximal in general
    keep='somatic*';     
    keepname='somatic';
    funcname='getKeepSomatic';
end
function [keep, keepname, funcname] = getKeepSomaticKCC2()
    keep='somatic_KCC2*';            
    keepname='somatic KCC2';
    funcname='getKeepSomaticKCC2';
end
function [keep, keepname, funcname] = getKeepSomaticNoKCC2()
    keep='somatic_[^K][^C]+[^2][^_].*';  
    keepname='somatic no KCC2';
    funcname='getKeepSomaticNoKCC2';
end

%% ANALYSIS FUNCTIONS
function [result, fig] = followDesired...
    (desiredFR,inhibition_strengths,excitation_strengths,firing_rates,xinhibit,titleName,units)
    if(nargin<6)
        titleName='';
        units='(Hz)';
    end
    if(strcmp(units(1),'(')~=1)
        units=strcat('(',units,')');
    end
    if(length(titleName)>1)
       titleName=strcat(titleName,'. '); 
    end
    
    for IN = 1:length(inhibition_strengths)
       spikelist = firing_rates{IN};
       [spikes, idx] = findClosestValue(desiredFR,spikelist,1);
       
       excite = excitation_strengths(idx);
       inhibit = str2double(inhibition_strengths{IN}(2:end));
       if(max(spikelist)<desiredFR)
           excite=200;              % inhibition is too strong, 
                                    % require more excitation 
                                    % to reach desired output rate
       end
       result{IN} = [excite inhibit];
       confidence(IN) = ((spikes-desiredFR));
    end
    exciteVals = cellfun(@(v) v(1), result(1,:));
    inhibitVals = cellfun(@(v) v(2), result(1,:));
    fig = figure();
    if(xinhibit==1)
        errorbar(inhibitVals,exciteVals,confidence,'-bo'); 
        plot(inhibitVals,exciteVals,'-bo'); 
        hold on;
        for res = 1:length(result)
           plot(result{res}(2),result{res}(1),'b*'); 
        end
        title(strcat(titleName,'Output rate: ',num2str(desiredFR),units));
        ylabel(strcat('Excitation ',units));
        xlabel(strcat('Inhibition ',units));
        try
            set(gca,'ylim',[0,exciteVals(end)]);
        catch
            % While setting the 'YLim' property of 'Axes':
            % Value must be a 1x2 vector of numeric type in which the second
            % element is larger than the first and may be Inf
            set(gca,'ylim',[0,1]);
        end
        try
            set(gca,'xlim',[0,inhibitVals(end)]);
        catch
            set(gca,'xlim',[0,1]);
        end
    else
        plot(exciteVals,inhibitVals,'-bo'); 
        hold on;
        for res = 1:length(result)
           plot(result{res}(1),result{res}(2),'b*'); 
        end
        title(strcat(titleName,'Output rate: ',num2str(desiredFR),units));
        xlabel(strcat('Excitation ',units));
        ylabel(strcat('Inhibition ',units));
        try
            set(gca,'xlim',[0,exciteVals(end)]);
        catch
            set(gca,'xlim',[0,1]);
        end
        try
            set(gca,'ylim',[0,inhibitVals(end)]);
        catch
            set(gca,'ylim',[0,1]);
        end
        
    end
end
function [value,idx] = findClosestValue(value,matrix,num)
    tmp = abs(matrix-value);
    if(num>0)
        idx = find(tmp==min(tmp),num);
    else
        idx = find(tmp==min(tmp));
    end
    
    value = matrix(idx);
end
function [Nindex1, Nindex2] = matchNindices(Nindex1,Nindex2)
    doSomething = 1;
    if(length(Nindex1)==length(Nindex2))
        doSomething = 0;
    end
    if(doSomething)
        if(length(Nindex1)>length(Nindex2))
            %shorten Nindex1
            Nindex1 = Nindex1(1:end-(length(Nindex1)-length(Nindex2)));
            firstN = Nindex1(1);
            lastN = Nindex1(end);
            newDif = firstN-lastN;
            Nindex1=(newDif-(1-Nindex1))./newDif;
        elseif(length(Nindex1)<length(Nindex2))
            %shorten Nindex2
            Nindex2 = Nindex2(1:end-(length(Nindex2)-length(Nindex1)));
            firstN = Nindex2(1);
            lastN = Nindex2(end);
            newDif = firstN-lastN;
            Nindex2=(newDif-(1-Nindex2))./newDif;
        end
    end
end

%% FITTING
function [params,stat,x,y] = sigmoid_fit(EX,WHY,fixed_params,initial_params,plot_flag)
    % fittings can change slightly depending on inclusion of fixed_params or
    % not ([])
    lastwarn(''); % reset warning state
    try
        [params,stat,x,y] = sigm_fit(EX,WHY,fixed_params,initial_params,plot_flag);
    catch Error
        disp(Error)
    end
     % check for error
    [warnmsg, msgid] = lastwarn;
    if ~strcmp(msgid,'')
      lastwarn(''); % reset warning state
      [params,stat,x,y] = sigm_fit(EX,WHY,[],initial_params,plot_flag);
      [warnmsg, msgid] = lastwarn;
      % check if there was still an error
      if ~strcmp(msgid,'')
          lastwarn(''); % reset warning state
          % try fitting with LESS parameters
          mask = WHY>0;
          [params,stat,x,y] = sigm_fit(EX(mask),WHY(mask),[],[],plot_flag);
      end
    end
%     mgompertz = @(p,x)(p(1)*exp(-exp(p(2)-p(3)*x)+p(4)*x));
%     [B,resnorm] = fminsearch(@(b) norm(WHY - mgompertz(b,EX)), [max(y) 1 1 1]);
%     ps = nlinfit(EX, WHY, gompertz, [max(y) 1 1 1]);
%     options=optimset('MaxFunEvals',10000,'MaxIter',5000);
%     [p,error]=lsqcurvefit(@gompertz,p0,EX,WHY,[],[],options);
end


%% QUANTIFY DIFFERENCES
function [Nindex, fig, subfigs] = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree,varargin)
    %% Difference BETWEEN 2 sets of parameter values
    % i.e. gives the difference between sets at parameter values.
    % For a given inhibition strength, how far is the midpoint of the 
    % IO curves (with and without chloride dynamics) from the baseline 
    % (no inhibition strength).
    % EX: vector of excitation strength
    % IN: vector in inhibition strength
    % INxFRFixed: vector of firing rate where chloride was fixed/static
    % INxFRFree: vector of firing rate where chloride was free/dynamic
    global toPlot toSave
    defaultOptions = struct('titlename','',...
                            'units','relative conductance',...
                            'plottype','semilogx',...
                            'plot',1,... % can explicitly suppress plotting
                            'ignore',[],...
                            'format','%.2f',...
                            'special',{''},...
                            'color','Blues',...
                            'synapses',0,...
                            'Nindex_fit', 0);
    options = processOptions(defaultOptions,varargin);
    if ( isa(options.units,'char'))
        options.units = options.units;
    end
    options.units = strcat('(',options.units, ')');
    INxEXFixed = 0;
    INxEXFree = 0;
    if (isa(EX,'cell'))
        INxEXFixed = EX{1};
        INxEXFree = EX{2};
        EX = INxEXFixed{1};
    end
    
    if(EX(1)==0)
        EX(1)=min(EX(2:end));
    end
    logEX = log10(EX);
    
    IN0 = IN(1);
    slopeBaseFixed = 1;
    x50BasedFixed = 0.5;
    special_points = [];
    subfigs=[];
    for k=1:length(INxFRFixed)
        if(options.plot>1)
            subfigs(k) = figure();
        end
        if (isa(INxEXFixed,'cell'))
            if(INxEXFixed{k}(1)==0)
                INxEXFixed{k}(1)=min(INxEXFixed{k}(2:end));
            end
            logEX = log10(INxEXFixed{k});
        end
        
        % Static/Fixed chloride
        fixed_params=[INxFRFixed{k}(1), INxFRFixed{k}(end), NaN, NaN];
        if(isnan(slopeBaseFixed))
            initial_params = [];
        else
            % set initial x50 to x050f (the base x50) and slope to slope0 (the base
            % slope), if this is not the first (base) run of course
            initial_params=[INxFRFixed{k}(1), INxFRFixed{k}(end), x50BasedFixed, slopeBaseFixed];
        end
        
        try
            [params,stat,x,y] = sigmoid_fit(logEX,INxFRFixed{k},fixed_params,initial_params,0);
            fitted_vals{k,1} = {x,y};
            % x50 value (X value where firing rate is half maximum)
            x50 = params(3);
            slope = params(4);
            if(options.plot>1)
                plots(1) = semilogx(10.^x,y,'b'); % points with fitted sigmoidal shape
                hold on;
                semilogx(INxEXFixed{k},INxFRFixed{k},'bo');
            end
            % from fitted values
            tmp = abs(x-x50);
            idx = find(tmp==min(tmp),1);
            if(options.plot>1)
                plot(10^x50,y(idx),'xr'); % red x on plot
            end
            if(k==1)
               % get base values (i.e. static chloride's lowest - hopefully 0 - inhibition)
               slopeBaseFixed=slope;
               xBaseFixed=x;
               yBaseFixed=y;
               x50BasedFixed=x50;
               idxBaseFixed=idx;
               EX0f=EX;
               INxFR0f = INxFRFixed{k};
            end
            % plot base
            if(options.plot>1)
                semilogx(10.^xBaseFixed,yBaseFixed,'k');
                semilogx(EX0f,INxFR0f,'ok');
                plot(10^x50BasedFixed,yBaseFixed(idxBaseFixed),'^k');
            end
        catch Error
            fitted_vals{k,1} = {nan, nan};
            disp(Error);
            disp(k);
            disp(logEX);
            disp(INxFRFixed{k});
            
        end
        
        
        % find firing rate (y) index of x50 value
        xq = logspace(log10(INxEXFixed{k}(1)), log10(INxEXFixed{k}(end)), length(INxEXFixed{k})*10);
        try
        yq = interp1(INxEXFixed{k},INxFRFixed{k},xq,'spline');
        catch Error
            xq = INxEXFixed{k};
            yq = INxFRFixed{k};
        end
        y50 = yq(end)/2;
        y_tmp = abs(yq-y50);
        y_idx = find(y_tmp==min(y_tmp),1);
        x50_ = xq(y_idx);
        
        if options.Nindex_fit
            pointfixed{k} = [10^x50 y(idx)];
        else
            pointfixed{k} = [x50_ yq(y_idx)];
        end
        %disp( pointfixed{k});
        
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Dynamic/Free chloride
        %%%%%%%%%%%%%%%%%%%%%%%%%
        fixed_params=[INxFRFree{k}(1), INxFRFree{k}(end), NaN, NaN];
        initial_params=[INxFRFree{k}(1), INxFRFree{k}(end), x50BasedFixed, slopeBaseFixed];
        
        if (isa(INxEXFree,'cell'))
            if(INxEXFree{k}(1)==0)
                INxEXFree{k}(1)=min(INxEXFree{k}(2:end));
            end
            logEX = log10(INxEXFree{k});
            
        end
        try
            [params,stat,x,y] = sigmoid_fit(logEX,INxFRFree{k},fixed_params,initial_params,0);
            fitted_vals{k,2} = {x,y};
            x50 = params(3);
            slope = params(4);
            if(options.plot>1)
                plots(2) = semilogx(10.^x,y,'g'); % points with fitted sigmoidal shape
                semilogx(INxEXFree{k},INxFRFree{k},'go');
            end
            tmp = abs(x-x50);
            idx = find(tmp==min(tmp),1);
            if(options.plot>1)
                plot(10^x50,y(idx),'xr'); % red x on plot
            end
            if(k==1)
               x0=x;
               y0=y;
               x050=x50;
               idx0=idx;
               EX0=EX;
               INxFR0 = INxFRFree{k};
            end
            % plot chloride base
            if(options.plot>1)
                plots(3) = semilogx(10.^x0,y0,'k');
                semilogx(EX0,INxFR0,'ok');
                plot(10^x050,y0(idx0),'*k');
                xlim('auto');
                % plot lines between base midpoint and latest inhibition
                line([10^x50BasedFixed,pointfixed{k}(1)],[yBaseFixed(idxBaseFixed),pointfixed{k}(2)]);
                line([pointfree{k}(1),pointfree{1}(1)],[pointfree{k}(2),pointfree{1}(2)]);
                %plotOrder = get(fig,'gca',plot);
                leg = legend(plots,{'Static Cl^-','Dynamic Cl^-','0 inhibition'});
                set(leg,'Location','NorthWest');
                inhTitle = IN(k);
                if(iscell(inhTitle))
                   inhTitle = inhTitle{1}; 
                end
                title(strcat(options.titlename, ': Inhibition: ',inhTitle(2:end), options.units));
                xlabel('Excitation');
                ylabel('Firing Rate (Hz)');
                set(gca,'xlim',[EX(1),EX(end)]);
                set(gca,'ylim',[0,max([max(y0),max(INxFRFixed{k}),max(INxFRFree{k})])]);
                if(and(toSave,toPlot))
                   saveFigs('folderName',{strcat(options.titlename, ': Inhibition: ',inhTitle(2:end), options.units)},subfigs(k));
                end
            end
        catch Error
            fitted_vals{k,2} = {nan, nan};
            disp(Error);
        end

        
        xq = logspace(log10(INxEXFree{k}(1)), log10(INxEXFree{k}(end)), length(INxEXFree{k})*10);
        try
            yq = interp1(INxEXFree{k},INxFRFree{k},xq,'spline');
        catch Error
            xq = INxEXFree{k};
            yq = INxFRFree{k};
        end
        y50 = yq(end)/2;
        y_tmp = abs(yq-y50);
        y_idx = find(y_tmp==min(y_tmp),1);
        x50_ = xq(y_idx);
        
        
        if options.Nindex_fit
            pointfree{k} = [10^x50 y(idx)];
        else
            pointfree{k} = [x50_ yq(y_idx)];
        end
        %disp( pointfree{k});
               
        
        if ismember(num2str(toNumber(IN{k})),options.special)
            try
                pointfixed{k}(2)
            catch Error
                disp(Error)
            end
            special_points = [special_points;[10^x50BasedFixed,pointfixed{k}(1)] [yBaseFixed(idxBaseFixed),pointfixed{k}(2)] [pointfree{k}(1),pointfree{1}(1)] [pointfree{k}(2),pointfree{1}(2)]];
        end
    end
    
    % calculate Normalisation index
    if options.Nindex_fit 
        for p=1:length(pointfixed)
           pointfixed{p} = log(pointfixed{p});
           pointfree{p} = log(pointfree{p}); 
        end
    end
    basePointFixed = pointfixed{1};
    basePointFree = pointfree{1};
    difWithIN0 = pdist([basePointFixed; basePointFree],'euclidean');
    disp(strcat(options.titlename,' difWithIN0= ',num2str(difWithIN0)));
    % check there was a fit
    for n=2:length(INxFRFixed)
        if(and(pointfixed{n}(1)>0,pointfree{n}(1)>0))
            fixedDist = pdist([basePointFixed; pointfixed{n}],'euclidean');
            freeDist = pdist([basePointFree; pointfree{n}],'euclidean');
            betweenPoints = pdist([pointfixed{n}; pointfree{n}],'euclidean');
            pfix = pointfixed{n};
            pfree = pointfree{n};
            pfree_2 = [pfree(1) pfix(2)];
            freeDist_2 = pdist([basePointFree; pfree_2],'euclidean');
            
            Nindex_orig = (fixedDist-freeDist)/fixedDist; % == 1-freeDist/fixedDist;
            Nindex_2 = 1-freeDist_2/fixedDist;
            Nindex(n) = Nindex_orig;
        else
            Nindex(n) = NaN;
        end
        if(options.plot>1)
            fig = figure(subfigs(n));
            text(100,100,strcat('NI:',num2str(Nindex(n)),' NI2:',num2str(Nindex_2)));
        end
    end
    
    if(options.plot>0)
        % plot fitted values in a nice single graph
        fig = figure();
        colors = cbrewer('seq', options.color, length(INxFRFixed));
        
        colors = [[0 0 0];colors(2:end,:)]; % start with black with 0 Inh
        color_index=1;
        special_lines = [];
        oranges = cbrewer('seq', 'YlOrRd', 2+length(options.special)); % minimum number is 3
        mid_Orange = oranges(round(length(oranges)/2),:);
        
        for n=1:length(INxFRFixed)
            if ismember(num2str(toNumber(IN{n})),options.special)
                color = oranges(length(special_lines)+2,:);
                lineWidth = 3;
            else
                color = colors(color_index,:);
                lineWidth = 2;
            end
            % plot dynamic
            
            plotsforlegend(n) = semilogx(10.^fitted_vals{n,2}{1},fitted_vals{n,2}{2},'-',...
                'Color',color,...
                'LineWidth',lineWidth,...
                'MarkerFaceColor', color);
            hold on;
            box off;

            % plot static
            ignoredplots(n)=semilogx(10.^fitted_vals{n,1}{1},fitted_vals{n,1}{2},'--',...
                'Color',color,...
                'LineWidth',lineWidth);


            semilogx(INxEXFree{n},INxFRFree{n},'x','Color',color);
            semilogx(INxEXFixed{n},INxFRFixed{n},'o','Color',color);

            color_index=color_index+1;
            if n==1
               example_dyn_static = [plotsforlegend(1).copy(), ignoredplots(1)];
            end

            if ismember(num2str(toNumber(IN{n})),options.special)
                special_lines = [special_lines;n];
            end

        end
        
        
        figure(fig.Number);
        for s_l=special_lines
           uistack(plotsforlegend(s_l),'top') 
           uistack(ignoredplots(s_l),'top')
        end
        
        uistack(plotsforlegend(1),'top') 
        uistack(ignoredplots(1),'top')
        
        num_special_points = size(special_points);
        for special_point_index=1:num_special_points(1)
            color = oranges(special_point_index+1,:);
            temp = special_points(special_point_index,:);
            x1 = temp(1);x2 = temp(2);
            y1 = temp(3);y2= temp(4);
            x3 = temp(5);x4 = temp(6);
            y3 = temp(7);y4 = temp(8);
            line([x1,x3],[y1,y3],                 'Color',color,'LineWidth',3,'MarkerSize',14);
            line([x1,x2],[y1,y2],'LineStyle','--','Color',color,'LineWidth',3,'MarkerSize',14)
            plot(x1,y1,'^','MarkerFaceColor',color,'Color',color,'MarkerEdgeColor','none','MarkerSize',14);
            plot(x2,y2,'x','MarkerFaceColor',color,'Color',color,'MarkerSize',14);
            plot(x3,y3,'x','MarkerFaceColor',color,'Color',color,'MarkerSize',14);
        end
        
        title(strcat(options.titlename));
        xlabel(strcat('Excitation ',' ',options.units));
        ylabel('Firing Rate (Hz)');
        YL = ylim;
        ylim([0,YL(2)]);
%         if options.synapses == 1
%             xlim([EX(1),50000]);
%         else
%             xlim([EX(1),1000]);
%         end
        % create new axis, and hide it
        ah=axes('position',get(gca,'position'));
        set(ah,'visible','off');
        leg_inh = legend(ah,example_dyn_static,{'Dynamic Cl^-','Static Cl^-'},'location','northeastoutside');
        legend boxoff
        numbersLegend = num2str(str2double(multiReplaceAll(IN,{'pt','d','c'},{'.','',''})),options.format);
        leg = legend(plotsforlegend,numbersLegend,'Location','EastOutside');
        
        t = title(leg,[{'Inhibition'} options.units]); % options.units
        legend boxoff
        
        % plot individual plots
        linestyles = {'--','-'};
        markerstyles = {'o','o'};
        INxFR = {INxFRFixed, INxFRFree};
        INxEX = {INxEXFixed, INxEXFree};
        
        fig_comb = figure();
        fig_saves = [fig_comb];
        % plot static and dynamic Cl i-o curves
        for kcc2=1:2
            fig_indv = figure();
            fig_saves = [fig_saves fig_indv];
            fig_nums = [fig_indv.Number fig_comb.Number];
            for f=1:length(fig_nums)
                figure(fig_nums(f));
                colors = cbrewer('seq', options.color, length(INxFRFixed));
                colors = [[0 0 0];colors(2:end,:)]; % start with black with 0 Inh
                color_index=1;
                special_lines = [];
                oranges = cbrewer('seq', 'YlOrRd', 2+length(options.special)); % minimum number is 3
                mid_Orange = oranges(round(length(oranges)/2),:);
                for n=1:length(INxFRFixed)
                    if ismember(num2str(toNumber(IN{n})),options.special)
                        color = oranges(length(special_lines)+2,:);
                        lineWidth = 3;
                    else
                        color = colors(color_index,:);
                        lineWidth = 2;
                    end
                    % plot either static or dynamic
    %                 temp_line(n) = semilogx(10.^fitted_vals{n,kcc2}{1},fitted_vals{n,kcc2}{2},linestyles{kcc2},...
    %                     'Color',color,...
    %                     'LineWidth',lineWidth);
    %                 hold on;
    %                 box off;
                    temp_line(n) = semilogx(INxEX{kcc2}{n},INxFR{kcc2}{n},...
                        strcat(linestyles{kcc2},markerstyles{kcc2}),...
                        'Color',color,...
                        'LineWidth',lineWidth,...
                        'MarkerFaceColor', color);
                    hold on;
                    box off;
                    color_index=color_index+1;

                    if ismember(num2str(toNumber(IN{n})),options.special)
                        special_lines = [special_lines;n];
                    end
                end
                if f==1 || kcc2==1
                    xlabel(strcat('Excitation',' ',options.units));
                    ylabel('Firing Rate (Hz)');
                    title(options.titlename);
                    YL = ylim;
                    ylim([0, YL(2)]);
                    numbersLegend = num2str(str2double(multiReplaceAll(IN,{'pt','d','c'},{'.','',''})),options.format);
                    leg = legend(temp_line,numbersLegend,'Location','EastOutside');    
                    t = title(leg,[{'Inhibition'} options.units]); % options.units
                    
                    % move special lines to top of z-index
                    for s_l=special_lines
                       uistack(temp_line(s_l),'top') 
                       uistack(temp_line(s_l),'top')
                    end

                    uistack(temp_line(1),'top') 
                    uistack(temp_line(1),'top')

                    num_special_points = size(special_points);
                    for special_point_index=1:num_special_points(1)
                        color = oranges(special_point_index+1,:);
                        temp = special_points(special_point_index,:);
                        x1 = temp(1);x2 = temp(2);
                        y1 = temp(3);y2= temp(4);
                        x3 = temp(5);x4 = temp(6);
                        y3 = temp(7);y4 = temp(8);
                        line([x1,x3],[y1,y3],                 'Color',color,'LineWidth',3,'MarkerSize',14);
                        line([x1,x2],[y1,y2],'LineStyle','--','Color',color,'LineWidth',3,'MarkerSize',14)
                        plot(x1,y1,'^','MarkerFaceColor',color,'Color',color,'MarkerEdgeColor','none','MarkerSize',14);
                        plot(x2,y2,'x','MarkerFaceColor',color,'Color',color,'MarkerSize',14);
                        plot(x3,y3,'x','MarkerFaceColor',color,'Color',color,'MarkerSize',14);
                    end

                end
            end
        end
        figure(fig_comb.Number)
        ah=axes('position',get(gca,'position'));
        set(ah,'visible','off');
        leg_inh = legend(ah,example_dyn_static,{'Dynamic Cl^-','Static Cl^-'},'FontSize', 22);
        legend boxoff
        
        if(and(toSave,toPlot))
           saveFigs('folderName',{strcat(options.titlename,'combined_fitted')},(fig));
           f_names = {'combined','static', 'dynamic'};
           for f=1:length(fig_saves)
               saveFigs('folderName',{strcat(options.titlename,f_names{f})},(fig_saves(f)));
           end
            
        end
    end
    
    if(~options.plot)
        % assign return argument that was not applicable
        fig = NaN;
        subfigs=NaN;
    end
end
function [Nindices, figs, fig] = CompareInnerNormalizationIndices...
    (Xmatrix,Ymatrix,rowentries,columnentries,figtitle,varargin)
    global toPlot toSave
    % set up default keyword arguments
    options = struct('orderby','column','rowtitle','','columntitle','');
    optionNames = fieldnames(options);
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('getCompare needs propertyName/propertyValue pairs')
    end
    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
       inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    
    if(or(strcmp(options.orderby,'column'),strcmp(options.orderby,'split')))
        Xmatrix=Xmatrix';
        Ymatrix=Ymatrix';
        fields=columnentries{1};
        legendentries=rowentries;
        legtitle=options.rowtitle;
        fieldtitle=options.columntitle;
    else
        fields=rowentries;
        legendentries=columnentries{1};
        legtitle=options.columntitle;
        fieldtitle=options.rowtitle;
    end
    for index = 1:length(fields)
         figtitlearg = strcat(figtitle,'(',fields{index},')');
         if(iscell(legendentries{1}))
             passlegend = legendentries{index};
         else
             passlegend = legendentries;
         end
         [Nindices{index},figs(index)] = ...
             InnerNormalizationIndex(Xmatrix(index,:),Ymatrix(index,:),...
             figtitlearg,passlegend,legtitle);
    end
     % plot normalisation indices against each other
     if(and(~isempty(Nindices),and(toPlot,~isnan(Nindices{1}))))
        fig = figure;
        hold on;

        for n = 1:length(Nindices)
            plotN=Nindices{n}(~isnan(Nindices{n}));
            %plot(legendentries{n}(1:length(plotN)),plotN,'-+');
            plot(1:length(plotN),plotN,'-+');
        end
        l = legend(fields);
        title(l,fieldtitle);
        set(gca,'XTick',linspace(1,length(legendentries),length(legendentries)));
        set(gca,'XTickLabel',legendentries);
        xlabel(legtitle);
        ylabel('"Normalisation" index');
        ytickformat('%.1f');
        figtitlearg = strcat(figtitle,'(',options.orderby,')');
        title(figtitlearg);
        if(toSave)
            saveFigs('folderName',{figtitlearg},(fig));
        end
     else
         fig = -1;
     end
end
function [Nindex, savefig] = InnerNormalizationIndex(EX,INxFR,figtitle,legendentries,legtitle)
    global toPlot toSave
    %% Difference WITHIN a set of parameter values
    % i.e. gives the difference between parameter values for a set.
    % To compare sets, may be able to compare returned Nindex from each
    % set.
    if(nargin<3)
       figtitle=''; 
    end
    if(nargin<4)
       legendentries={}; 
    end
    if(nargin<5)
        legtitle='';
    end
    if(toPlot)
        savefig=figure();
    else
        savefig=-1;
    end
    for k=1:length(INxFR)
        if(isempty(INxFR{k}))
           continue; 
        end
        % convert to log-friendly value
        if(EX{k}(1)==0)
            EX{k}(1)=1;
        end
        % get log values
        logEX = log10(EX{k});
        maxlogEX = max(logEX);
        % set sigmoidal parameters (min, max, mid, slope)
        if(length(logEX)>2)
            initial_params=[min(INxFR{k}), max(INxFR{k}) , NaN , NaN];
        else
            initial_params=[NaN, NaN , NaN , NaN];
        end
        % fit to sigmoid
        try
            [params,stat,x,y] = sigm_fit(logEX,INxFR{k},[],initial_params,0);
        catch
            [params,stat,x,y] = sigm_fit(logEX,INxFR{k},[],[],0);
        end
        x50 = params(3);
        slope = params(4);
        % plot result
        if(toPlot)
            % fitted sigmoid
            plots(k)=semilogx(10.^x,y);
            c = get(plots(k),'Color');
            hold on;
            % actual points
            temp = semilogx(EX{k},INxFR{k},'-.o');
            set(temp,'Color',c);
        end
        % check if x50 is to the far right of the graph (happens on extreme
        % values of non-sigmoidal fits)
        if(x50>maxlogEX)
           %10^x50 == Inf
           x50=logEX(end);
           idx = length(y);
        else
            % determine each point's distance away from middle
            tmp = abs(x-x50);
            % y-index of x50
            idx = find(tmp==min(tmp),1);
        end
        % plot x50 point
        if(toPlot)
            plot(10^x50,y(idx),'xr');
        end
        % store x50
        EX50(k) = x50;
        % and corresponding (fitted) firing rate
        FR50(k) = y(idx);
        % as a point
        point{k} = [EX50(k) FR50(k)];
        fprintf('EX50:%f FR50:%2.2f\n',EX50(k),FR50(k));
    end
    % make the plot prettier
    if(toPlot)
        title(figtitle);
        leg = legend(plots,legendentries);
        set(leg,'Location','NorthWest');
        title(leg,legtitle);
        xlabel('Excitation');
        ylabel('Firing Rate');
    end
    
    % remove infinities
    EX50=EX50(~isinf(EX50));
    % find max x50
    maxEX50 = max(EX50);
    % get index of maximum
    match = find(EX50==maxEX50);
%     if(match<length(EX50))
%        for i = match+1:length(EX50)
%           EX50(i) = EX(end);
%           point{i} = [EX50(i) FR50(i)];
%        end
%        match=length(EX50);
%     end
    % get corresponding index for firing rate
    FRmaxEX50 = FR50(match);
    % get corresponding index for point
    maxPoint = point{match};
    % set base as first point
    basePoint = point{1};
    % greatest distance is between base and max points
    maxDist = pdist([basePoint; maxPoint],'euclidean');
    for n=1:match   % ignore lesser EX50 points that are above match (usually due to flat line)
        Nindex(n) = (maxDist-pdist([basePoint; point{n}],'euclidean'))/maxDist;
    end
    if(and(toSave,toPlot))
       saveFigs('folderName',{figtitle},(savefig));
    end
end

%% RETRIEVE FROM DATA
function [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2...
        (foldername,extraCondition,prox,dist,inter,plot,save,plottype,varargin)
    options = struct('excludeSplit',{''},...
                     'excludeDegree',{''});
    options = processOptions(options, varargin);
    if(nargin<8)
        plottype='normal';
    end
    if(prox)
        keep=strcat('proximal.*',extraCondition);  % proximal in general
        [proximal,proximal_KCC2] = compareKCC2(foldername,keep,...
                                    'excludeSplit',options.excludeSplit,...
                                    'excludeDegree',options.excludeDegree); 
    else
        proximal = [];
        proximal_KCC2=[];
    end
    if(dist)
        keep=strcat('distal.*',extraCondition);
        [distal, distal_KCC2] = compareKCC2(foldername,keep,...
                                    'excludeSplit',options.excludeSplit,...
                                    'excludeDegree',options.excludeDegree);
    else
        distal = [];
        distal_KCC2=[];
    end
    if(inter)
        keep=strcat('internrn.*',extraCondition);
        [internrn, internrn_KCC2] = compareKCC2(foldername,keep,...
                                    'excludeSplit',options.excludeSplit,...
                                    'excludeDegree',options.excludeDegree);
    else
        internrn=[];
        internrn_KCC2=[];
    end
    
    if(or(plot,save))
        figs = [];
        printNames = {};
        if(prox)
            pFig = plotStruct(proximal,'Proximal (no chloride dynamics)','Excitation','Spikes','plottype',plottype);
            pkFig = plotStruct(proximal_KCC2,'Proximal (with chloride dynamics)','Excitation','Spikes','plottype',plottype);
            figs = [pFig pkFig];
            printNames={'Proximal no KCC2','Proximal KCC2'};
        end
        if(dist)
            dFig = plotStruct(distal,'Distal (no chloride dynamics)','Excitation','Spikes','plottype',plottype);
            dkFig = plotStruct(distal_KCC2,'Distal (with chloride dynamics)','Excitation','Spikes','plottype',plottype);
            figs = [figs dFig dkFig];
            combineFigures([dFig dkFig],'titlename','Distal',...
                            'plottype','logx',...
                            'matchcolors',1,...
                            'uniquespecifiers',0,...
                            'lines',1,...
                            'legend','once');
            printNames = {printNames{1:end},'Distal no KCC2','Distal KCC2','plottype','normal'};
        end
        if(inter)
            iFig = plotStruct(internrn,'Interneuron (no chloride dynamics)','Excitation','Spikes','plottype',plottype);
            ikFig = plotStruct(internrn_KCC2,'Interneuron (with chloride dynamics)','Excitation','Spikes','plottype',plottype);
            figs = [figs iFig ikFig];
            printNames = {printNames{1:end},'Interneuron no KCC2','Interneuron KCC2'};
        end
        if(save)
            printNames=strcat(printNames,extraCondition);
            saveFigs(foldername,printNames,figs);
        end
    end
    
end
function [noKCC2, withKCC2] = compareKCC2(foldername,keep,varargin)
    options = struct('excludeSplit',{''},...
                     'excludeDegree',{''});
    options = processOptions(options, varargin);
    output = getCompare(foldername,keep,'KCC2','InPt[.\d]+[(_]',...
        'trimBeginDegree',4,...
        'trimEndDegree',1,...
        'excludeSplit',options.excludeSplit,...
        'excludeDegree',options.excludeDegree);
    noKCC2 = output.sno;
    withKCC2 = output.sKCC2;
end

%% GET DATA AND PLACE IN STRUCT
function [resultContainer,colNames] = getCompare(foldername,keep,splitReg,degree,varargin)
    %# define defaults at the beginning of the code so that you do not need to
    %# scroll way down in case you want to change something or if the help is
    %# incomplete
    options = struct('split','match',...
                     'degree','match',...
                     'excludeSplit',{''},...
                     'excludeDegree',{''},...
                     'trimBeginSplit',0,...
                     'trimEndSplit',0,...
                     'trimBeginDegree',0,...
                     'trimEndDegree',0);
    options = processOptions(options, varargin);

    files = dir(strcat(foldername,'/*.dat')); %must be exact

    [matrixDimen,colNames,formatSpec,headerLines] = getMatrixDimen(strcat(foldername,'/',files(1).name));
    columns = str2double(matrixDimen{2});
    disp('Headers');
    disp(colNames);
    %# begin reading of each file
    for i=1:length(files)
        %# keep only subset
        strRegEx=files(i).name;
        if(isempty(regexp(strRegEx,keep, 'once')))
            continue;
        end
        
        %# read data
        fileID = fopen(strcat(foldername,'/',strRegEx),'r');
        indata = textscan(fileID,formatSpec,'HeaderLines',headerLines);
        fclose(fileID);
        
        %# place data in container struct 
        for col = 1:columns
            container.(colNames{col})= indata{col};
        end
        container.Spikes_at_Axon = container.Spikes;
%         container.Spikes = container.Spikes_at_Soma;
%         [argvalue, argmax] = max(container.Spikes);
%         if (argvalue <= 1)
%            argmax=length(container.Spikes);
%         end
%         container.Spikes = container.Spikes(1:argmax);
%         container.Excitation = container.Excitation(1:argmax);
        
        %# split files into categories
        %#  e.g. number of synapses
        disp(strRegEx);
        intermediateSplitResult = regexp(strRegEx,splitReg, 'once',options.split);
        intermediateSplitResult = intermediateSplitResult(1+options.trimBeginSplit:end-options.trimEndSplit);
        if(ismember(intermediateSplitResult,options.excludeSplit))
            % skip this one
            continue; 
        end
        splitResult = sanitizeResult(intermediateSplitResult,'start','s');
        if(strcmp(splitResult,'s'))
            splitResult = 'sno';
        end
        
        %# each category should have a 'degree' of another measure
        %#  e.g. inhibition strength
        if(strcmp(options.degree,'split'))
           degreeString = regexp(strRegEx,degree,options.degree);
           degreeTake = degreeString{end}; 
           degreeNumbers = regexp(degreeString{end},degree, 'match');
           degreeResult = strcat('d',multiReplace(sanitizeNum(degreeNumbers(1)),{'[.]' '_' '-'},{'pt' '' 'm'}));
        else
            intermediateDegreeResult = regexp(strRegEx,degree, 'once',options.degree);
            intermediateDegreeResult = intermediateDegreeResult(1+options.trimBeginDegree:end-options.trimEndDegree);
            if(ismember(num2str(str2double(intermediateDegreeResult)),options.excludeDegree))
                % skip this one
                continue; 
            end
            degreeResult = sanitizeResult(intermediateDegreeResult,'start','d');
        end
        
        disp(strcat('splitResult:',splitResult,' degreeResult:',degreeResult));
        resultContainer.(splitResult).(degreeResult) = container;

    end
    disp('Splits');
    disp(resultContainer);
end
%% TRANSPOSE
function new_output = transposeSynapticNumberToStrength(output)
    split_fields = fields(output);
    a = split_fields{1};
    c = strfind(a,'c');
    exc = a(1:c);
    % find how many elements contain the first excitation synapse
    % number value (i.e. excitation group size)
    exc_size = sum((strfind([split_fields{:}],exc)>0));
    % calculate inhibitory group size
    inh_size = length(split_fields)/exc_size;
    new_output = struct;
    for s = 1:length(split_fields)
        % find 'c' - middle of numbers
        c = strfind(split_fields{s},'c');
        % inhibitory number of synapses
        new_degree = split_fields{s}(c:end);
        % excitatory number of synapses
        excitation = toNumber((split_fields{s}(2:c-1)));

        disp(strcat('inh:',new_degree,'exc:',num2str(excitation)));
        splits = fields(output.(split_fields{s}));
        for split_index = 1:length(splits)
            split=splits{split_index};
            % number of spikes
            tmp_struct = output.(split_fields{s}).(split);
            for degree_num_index = 1:length(tmp_struct.Excitation)
                degree_num = tmp_struct.Excitation(degree_num_index);
                new_split = strcat('e',num2str(degree_num),'i',split(2:end));
                spikes = tmp_struct.Spikes(tmp_struct.Excitation == degree_num);
                other_fields = fields(tmp_struct);
                % add to matrix
                try
                    new_output.(new_split).(new_degree).Excitation = [new_output.(new_split).(new_degree).Excitation, excitation];
                    new_output.(new_split).(new_degree).Spikes = [new_output.(new_split).(new_degree).Spikes, spikes];
                    for other_field_index = 1:length(other_fields)
                        other_field = other_fields{other_field_index};
                        if or(strcmp(other_field,'Excitation'), strcmp(other_field,'Spikes'))
                            continue
                        else
                            new_output.(new_split).(new_degree).(other_field) = ...
                                [new_output.(new_split).(new_degree).(other_field), ...
                                tmp_struct.(other_field)(tmp_struct.Excitation == degree_num)];
                        end
                    end
                catch
                    % first encounter of new_split > new_degree
                    new_output.(new_split).(new_degree).Excitation(1) = excitation;
                    new_output.(new_split).(new_degree).Spikes(1) = spikes;

                    for other_field_index = 1:length(other_fields)
                        other_field = other_fields{other_field_index};
                        if or(strcmp(other_field,'Excitation'), strcmp(other_field,'Spikes'))
                            continue
                        else
                            new_output.(new_split).(new_degree).(other_field)(1) = ...
                                tmp_struct.(other_field)(tmp_struct.Excitation == degree_num);
                        end
                    end
                end
            end
        end
    end
    % everything is in the struct, but may be unordered.
    new_fields = fields(new_output);
    for new_fields_index = 1:length(new_fields)
        new_field = new_fields{new_fields_index};
        new_sub_fields = fields(new_output.(new_field));
        for new_sub_fields_index = 1:length(new_sub_fields)
            new_sub_field = new_sub_fields{new_sub_fields_index};
            tmp_new_struct = new_output.(new_field).(new_sub_field);
            recording_fields = fields(tmp_new_struct);
            excitation = tmp_new_struct.Excitation;
            for recording_field_index = 1:length(recording_fields)
                recording_field = recording_fields{recording_field_index};
                if strcmp(recording_field,'Excitation')
                    continue
                else
                    % sort everything according to Excitation (a
                    % persistently unordered excitation vector is used to
                    % reorder consistently)
                    [new_output.(new_field).(new_sub_field).Excitation,...
                        new_output.(new_field).(new_sub_field).(recording_field)]...
                        = multiSort(excitation,tmp_new_struct.(recording_field));
                end
            end
        end
    end
    
end
%% HELPER FUNCTIONS
function options = processOptions(options,inVarargin)
    %# read the acceptable names
    optionNames = fieldnames(options);
    
    %# count arguments
    nArgs = length(inVarargin);
    if round(nArgs/2)~=nArgs/2
       error('need propertyName/propertyValue pairs')
    end

    for pair = reshape(inVarargin,2,[]) %# pair is {propName;propValue}
       inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          %# overwrite options. If you want you can test for the right class here
          %# Also, if you find out that there is an option you keep getting wrong,
          %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
end
function [a,b] = multiSort(a,b)
    c = [a;b]';
    d = sortrows(c)';
    a=d(1,:);
    b=d(2,:);
end
function result = sanitizeResult(string,varargin)
    %% Convert name to a santized version MATLAB allows as dictionary index
    options = struct('start','',...
        'toReplace',{{'[.]' '_' '-' '[' ']' '(' ')'}},...
        'replacements',{{'pt' '' 'm' '' '' '' ''}});
    %# read the acceptable names
    optionNames = fieldnames(options);
    
    %# count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('plotAndSave needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
       inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          %# overwrite options. If you want you can test for the right class here
          %# Also, if you find out that there is an option you keep getting wrong,
          %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    string = multiReplace(string,{','},{'c'});  %replace this before next step as MATLAB ignores commas
    result = strcat(options.start,multiReplace(sanitizeNum(string),options.toReplace,options.replacements));
end
function sanNumber = sanitizeNum(number)
    %% Remove leading/trailing zeroes
    %if not a number return back as it was
    sanNumber = num2str(str2double(number));
    if(strcmp(sanNumber,'NaN'))
       sanNumber  =  number;
    end
end
function string = multiReplace(string,toReplace,replacements)
    %% more powerful version of regexprep
    %   the toReplace cell is matched with the replacements cell so that
    %   the 1st index of of toReplace is replaced by the 1st index of
    %   replacements, 2nd -> 2nd, 3rd -> 3rd.
    %   an uneven match will result in an error.
    %   example:
    %       multiReplace('the string to replace',{' ','to','replace'},...
    %                                            {'_','is','replaced'})
    %       ans = the_string_is_replaced
    for tr = 1:length(toReplace)
        string  = regexprep(string,toReplace(tr),replacements(tr));
    end
end
function newString = multiReplaceAll(string,toReplace,replacements)
    %% more powerful version of multiReplace
    %   all replacements must be made for a successful replacement
    newString = string;
    oldString = string;
    for tr = 1:length(toReplace)
        newString  = regexprep(oldString,toReplace(tr),replacements(tr));
        oldString = newString;
    end
end
function [matrixDimen,colNames,formatSpec,headerLines] = getMatrixDimen(fileStr)
    %% determine matrix dimensions and column names (tab delimited) from a file
    fileID = fopen(fileStr,'r');
    formatSpec='%s';
    indata = textscan(fileID,formatSpec,'Delimiter','\t');
    fclose(fileID);
    outPut = indata{1};
    headerLines=0;
    formatSpec='';
    for i=1:length(outPut)
       outP = outPut{i};
       if(~isempty(str2num(outP))) %#ok<ST2NM>
           s = colNames{i-1};
           matrixDimen = regexp(s,'[\d]*','match');
           if(length(matrixDimen)==2)
               % matrix dimensions on same line as headers
               headerLines=1;
           else
               % matrix dimensions on next line
               matrixDimen = regexp(outP,'[\d]*','match');
               headerLines=2;
           end
           colNames{i-1} = regexprep(s,'\d','');
           break;
       else
           colNames{i} = regexprep(strtrim(outP),'\s','_');
           formatSpec = strcat(formatSpec,'%f');
       end
    end
end
function result = toNumber(string,cutEnd)
    %% converts previosly sanitized number (see sanitizaResult) back to a number
    if(nargin==1)
       cutEnd = 0; 
    end
    string = multiReplace(string,{'d','s','c','pt','m'},{'','',',','.','-'});
    if(iscell(string))
        for s = 1:length(string)
            result(s) = str2double(string{s}(1:end-cutEnd));
        end
    else
        result = str2double(string(1:end-cutEnd));
    end
    
end
function [reordered,sanitized] = reorderByNum(stringNumbers,cutEnd)
    if(nargin==1)
       cutEnd = 0; 
    end
       order = 1:length(stringNumbers);
       % Convert string to number
       order = toNumber(stringNumbers,cutEnd);
       if isnan(order(1))
           order = zeros(length(stringNumbers),...
                         length(split('c',stringNumbers{1})));
           for n=1:length(stringNumbers)
               order(n,:) = toNumber(split('c',stringNumbers{n}),cutEnd);
           end
           % Sort and get new indices
           [neworder,newindices] = sort(order(:,1));
       else
           % Sort and get new indices
           [neworder,newindices] = sort(order);
       end
       % Reorder original
       reordered = stringNumbers(newindices);
       % convert numbers to string for returning sanitized cell
       for no = 1:length(neworder)
            sanitized{no} = num2str(neworder(no));
       end
end
function fitted_y = lineOfBestFit(x,y,poly)
    if(nargin==2)
        poly=1;
    end
    fitted_y = polyval(polyfit(x,y,poly),x);
end
%% PLOT
function colors = getColorPalette(varargin)
    % Matlab built in colour palettes here: 
    %   https://www.mathworks.com/help/matlab/ref/colormap.html
    % Default color palette
    % colors = get(groot,'DefaultAxesColorOrder');
    
    defaultOptions = struct('Palette','all',...% Great data viz palette
                            'Black',false,...
                            'N', 10); 
    options = processOptions(defaultOptions,varargin);
    try
        colors = palette(options.Palette);
    catch
        colors = cbrewer('seq', options.Palette, options.N);
    end
 
    if options.Black == 1
        colors = [[0 0 0];colors(2:end,:)]; % start with black
    end
end
function fig = plotStruct(struct1, name,fieldX,fieldY,varargin)
    %% plotStruct
    % struct1, name, fieldX, fieldY, varargin

    defaultOptions = struct('colorstring',getColorPalette(),...
                            'specifier','o+*.xsd^v><ph',...
                            'plottype','semilogx');
    options = processOptions(defaultOptions,varargin);
    
    
   %   plot graphs
   fig = figure();

   plotMat = zeros(length(fields(struct1)),1);
   
   e=1;
   ce=1;
   se=1;
   structFields = fields(struct1);
   for f = 1:length(structFields)
       field = structFields{f};
       if(strcmp(options.plottype,'semilogx'))
          if(struct1.(field).(fieldX)(1)==0)
              struct1.(field).(fieldX)(1)=1;
          end
          plotMat(e) = semilogx(struct1.(field).(fieldX),struct1.(field).(fieldY),strcat('-',options.specifier(se)),'Color',options.colorstring(ce,:));
          hold('on');
       else
          plotMat(e) = plot(struct1.(field).(fieldX),struct1.(field).(fieldY),strcat('-',options.specifier(se)),'Color',options.colorstring(ce,:));
          hold('on');
       end
       
       e=e+1;
       ce=ce+1;
       se=se+1;
       if(ce>length(options.colorstring))
           ce=1;
       end
       if(se>length(options.specifier))
           se=1;
       end
   end

   %order = names;


   % Set plot properties including title and legend
   try
       figure(fig);
       title(name);
       xlabel(fieldX);
       ylabel(fieldY);
       %% reorder plot and legend so 0 is first (happens when values (0-1) exist)
       numbersLegend = regexprep(regexprep(structFields,'d',''),'pt','.');
       order = 1:length(numbersLegend);
       if(strcmp(numbersLegend{1},'0')~=1)
          tempMatrix = strcmp(numbersLegend,'0');
          for i = 1:length(tempMatrix)
             if(tempMatrix(i)==1)
                 index=i;
                 order = [index 1:index-1 index+1:length(tempMatrix)];
                 break;
             end
          end
       end
       
       
       %% set legend
       leg = legend(plotMat(order),numbersLegend{order});
       try
           v = get(leg,'title');
           set(v,'string','Inhibition');
       catch ME
           disp(ME.message)
       end
       set(leg,'Location','EastOutside');


   catch ME 
       disp('legend error');
       disp(ME.message);
   end
    
end
function fig = formatFigure(fig,varargin)
    defaultOptions = struct('PaperSize',[8.27,11.69],...
                            'ImageScale',[1,0.5]); % value between 0 and 1 (0-100%)
    options = processOptions(defaultOptions,varargin);
    % make a figure pretty
    if length(isgraphics(fig))>1
        for f = 1:length(fig)
            fig(f) = formatFigure(fig(f));
        end
    else
        % get parent axes
        axes1 = get(fig,'CurrentAxes');
        % Create multiple lines using matrix input to plot
        set(findall(fig, 'Type', 'Line'),'LineWidth',2);
        % Set the remaining axes properties
        set(axes1,'FontSize',14);
        % Format legend
        legend1 = legend(axes1);
        set(legend1,'Location','northwest','FontSize',14);
        set(fig, 'Units', 'Inches', 'Position', [[0,0], options.ImageScale.*options.PaperSize],...
            'PaperUnits', 'Inches', 'PaperSize', options.PaperSize)
    end
end
function [fig1,fig2] = adjustLims(fig1,fig2)
    %% Set figures to have the same x and y limits
    %   TODO: allow variable number figures
    Ax = get(fig1,'CurrentAxes');
    kAx = get(fig2,'CurrentAxes');
    y1 = get(Ax,'ylim');
    y2 = get(kAx,'ylim');

    if(y2(2)>y1(2))
      y1(2) = y2(2);
    else
      y2(2) = y1(2);
    end

    set(Ax,'ylim',y1);
    set(kAx,'ylim',y2);
end
function newfig = combineFigures(h,varargin)
    %% Place 2 figures into a combined figure such that data can be better compared
    % h should be a Vector of figs such that h(1) = figure(); h(2) =
    % figure();
    % use 'plottype','logx' to use semilogx plot type
    defaultOptions = struct('titlename','',...
                            'plottype','logx',...
                            'matchcolors',0,...
                            'uniquespecifiers',1,...
                            'lines',1,...
                            'legend','full');   % full or once
    options = processOptions(defaultOptions,varargin);
    
    newfig = figure();
    box off;
    handleLines = findobj(h,'type','line');
    handleAxes = findobj(h,'type','axes'); % can get axis scale type (log/linear) from these objects
    legends = findobj(h,'type','legend');
    if(isempty(legends))
        legLen = length(handleLines);
    else
        % need to flip as they are inverse of handles
        l = flipud(legends(1)); % assuming legends all the same...
        legLen = length(l.String);
    end
    hLen = length(h);
    lLen = length(handleLines);
    % if not every line was in the legend, adjust number of legend entries
    while(and(mod(lLen,legLen)~=0,lLen>=legLen))
        legLen=legLen+1;
    end
    legIndices=[];
    legIndicesIndex = 1;
    hIndex = 1;
    baseTitle = handleAxes(hIndex).Title.String;
    lstringIndex=1;
    
    % get correct display names for each line
    for i=1:lLen
        if(~isempty(handleLines(i).DisplayName))
            legIndices(legIndicesIndex) = i;
            % name plus figure origin
            legFull{legIndicesIndex} = strcat(baseTitle,' - ',handleLines(i).DisplayName);
            if(i<=lLen/hLen)
                % just name (original legend entry)
                legOnce{legIndicesIndex} = handleLines(i).DisplayName;
            end
            if(and(i/lLen==hIndex/hLen,hIndex<hLen))
                % next source figure/axis
                hIndex=hIndex+1;
                baseTitle = handleAxes(hIndex).Title.String;
                lstringIndex=0;
            end
            legIndicesIndex=legIndicesIndex+1;
            lstringIndex=lstringIndex+1;
        end
    end
    
    % set color palette  and specifiers
    colors = getColorPalette();
    specifier = 'osd+x*.^v><ph';
    c=1;
    if(options.lines)
        s=0;
    else
        s=3; 
    end
    
    plotsIndex=1;
    % iterate over source figure lines (in groups of figures' lines)
    % such that each i is the index of the first line from a source figure
    for i = (lLen/hLen)/legLen:(lLen/hLen):lLen
        if(options.uniquespecifiers)
            if(options.lines) 
                s=1;
            else
                s=3; 
            end
        else
            s=s+1;
        end
        
        if(options.matchcolors)
            % have the colors start from 1 for each source figure
            c=1;
        end
        for legIndex = 0:legLen-1
            % go through each line of one of the source figures
            index = i+legIndex*((lLen/hLen)/legLen);
            index = max(1,round(index));
            if(options.lines)
                plotting_aesthetics = strcat('-',specifier(s));
            else
                plotting_aesthetics = specifier(s);
            end
            
            if(strcmp(options.plottype,'logx'))
                tmp = semilogx(handleLines(index).XData, handleLines(index).YData,...
                    plotting_aesthetics,'Color',colors(c,:));
            else
                tmp = plot(handleLines(index).XData, handleLines(index).YData,...
                    plotting_aesthetics,'Color',colors(c,:));
            end
            if(plotsIndex<=length(legIndices) && index==legIndices(plotsIndex))
                % note if this plot should be mentioned in the legend
                plots(plotsIndex) = tmp;
                plotsIndex=plotsIndex+1;
                if(options.matchcolors)
                    % increase colour per (legend) line
                    c=c+1;
                    if(c>length(colors))
                        c=1;
                    end
                end
            end
            % important to have hold after plotting of logx was chosen
            hold on;
            
            if(options.uniquespecifiers)
                s=s+1;
                if(s>length(specifier))
                    s=1;
                end
            end
        end
        if(~options.matchcolors)
            % increase colour per source figure (each source figure will have different
            % colour)
            c=c+1;
            if(c>length(colors))
                c=1;
            end
        end
    end
    
    % add legend
    if(strcmp(options.legend,'full'))
        % legend includes which source figure each line came from
        try
            legend(plots, legFull);
        catch
            try
                legend(legFull);
            catch
                legend();
            end
        end
    else
        % legend is each corresponding line from all the source figures
        try
            legend(plots, legOnce);
        catch
            try
                legend(legOnce);
            catch
                legend();
            end
        end
    end
    xlabel(handleAxes(1).XLabel.String);
    ylabel(handleAxes(1).YLabel.String);
    title(options.titlename);
    box off;
end
%% SAVE
function saveFigs(foldername,saveNames,figs,quick,varargin)
    global toSave quickPlot folderName overWrite
    if(toSave~=1)
       disp('ERROR: saveFigs called but toSave not 1');
       return
    end
    if(nargin<4)
      % if quick argument not passed, set to global value
      quick = quickPlot; 
    end
    if(strcmp(foldername,'folderName'))
        % check global variable if requested
        foldername=folderName;
    end
    old = cd(foldername);
    saveNames = strcat('graph_',saveNames);
    posX=0;
    posY=0;
    index=1;
    diff=3;
    if(quick)
       saveFormats = {'-dpng'};
    else
       saveFormats = {'-depsc','-dpng'};
    end
    figN = 1;
    for fig = figs
       for sf=1:length(saveFormats)
           saveFormat=saveFormats{sf};
           disp(fig);
    %                pnFig=fig;
    %                while(pnFig>length(printNames))
    %                    pnFig=pnFig-6;
    %                end
           savename = multiReplace(saveNames{figN},{':',' I',' '},{'',',I','_'});
           if ~overWrite
               savename_index=1;
               while exist(strcat(savename,'.png'), 'file') == 2
                   savename = strcat(savename,'_',num2str(savename_index));
                   savename_index = savename_index+1;
               end
           end
           % save figure
           print(fig,savename,saveFormat,'-r0');
           disp(strcat('saved ', savename,' ',saveFormat));   
       end
       figN=figN+1;
       set(fig, 'Position', [posX posY 1920/3 1080/2]);
       %disp(strcat('posX',num2str(posX)));
       %disp(strcat('posY',num2str(posY)));
       if(index==3)
           posX = 0;
           diff=diff*2;
       else
           posX=1920*(diff-index)/3;
       end

       if(index==3)
           posY = 1080/2;
       end

       set(fig,'PaperPositionMode','auto')
       index = index + 1;

    end

    cd(old);
end

%% COMMON ANALSYSES
function result = SplitOnNameDegreeOfInhibition(foldername, keepFunctions, varargin)
    defaultOptions = struct('NindexTypes','both');
    options = processOptions(defaultOptions,varargin);
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        [split_fields, split_fields_num] = reorderByNum(fields(output));

        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
        end
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        NindexIndex = 1; % confusing name much...
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'split')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','split','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'degree')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','degree','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end
    end
    result = 1; 
end

function result = CompareKeepsSplitOnNameDegreeOfInhibition(foldername, keepFunctions, varargin)
    defaultOptions = struct('NindexTypes','both');
    options = processOptions(defaultOptions,varargin);
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(funcnames{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));

        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
        end
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        topoutput.(funcnames{k}).excitation = excitation;
        topoutput.(funcnames{k}).spikes = spikes;
        NindexIndex = 1; % confusing name much...
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'split')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','split','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'degree')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','degree','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end  
    end
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(funcnames{k}).excitation{split_index,:},topoutput.(funcnames{k}).spikes{split_index,:})
        end
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    specifiers = ' -o+*.xsd^v><ph';
    if(length(funcnames)>length(specifiers))
       disp(strcat('ERROR in CompareKeepsSplitOnNameDegreeOfInhibition: ','length(funcnames)>length(specifiers)'));
       pause;
    end
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(funcnames{1}).excitation{split_index,:},topoutput.(funcnames{1}).spikes{split_index,:},'-','Color',colors(split_index,:));
        for funcIndex = 2:length(funcnames)
            plot(topoutput.(funcnames{funcIndex}).excitation{split_index,:},topoutput.(funcnames{funcIndex}).spikes{split_index,:},strcat('-',specifiers(funcIndex)),'Color',colors(split_index,:));
        end
    end
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    titlename = strcat(funcnames{1},':-');
    for funcIndex = 2:length(funcnames)
        titlename=strcat(titlename,' ',funcnames{funcIndex},strcat(':-',specifiers(funcIndex)));
    end
    title(titlename);
    saveFigs('folderName',{'keep functions compared'},f3);
    result=1;
end

function result = splitOnKCC2DegreeOfName(foldername, keepFunctions,varargin)
    defaultOptions = struct('NindexTypes','both');
    options = processOptions(defaultOptions,varargin);
    
    [split,splitName,trimBeginSplit,trimEndSplit]=getSplitOnKCC2();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfNameInPt(); 
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields_sub = multiReplaceAll(degree_fields_sub,{'d','c'},{'E:',' I:'});
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
        end
        split_fields = multiReplaceAll(split_fields,{'s'},{' '});
        NindexIndex = 1; % confusing name much...
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'split')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','split','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end
        if(or(strcmp(options.NindexTypes,'both'),strcmp(options.NindexTypes,'degree')))
            [Nindices{NindexIndex}, figs{NindexIndex}, fig{NindexIndex}] = CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                    degree_fields,keepname,'orderby','degree','rowtitle',...
                    splitName,'columntitle',degreeName);
            combinedfig(NindexIndex) = combineFigures(figs{NindexIndex});
        end
    end
    result = 1;
end

function [result,fig,subfigs] = NormalisationIndexAnalysisBase(data,data_KCC2,EX,IN,varargin)
        global toSave toPlot
        defaultOptions = struct('titlename','',...
                                'format','%.2f',...
                                'special',{''},...
                                'units','',...
                                'color','Blues');
        options = processOptions(defaultOptions,varargin);
        x_index=1;
        [IN,IN_num] = reorderByNum(fields(data));
        for x=1:length(IN)
            INxEXFixed{x} = data.(IN{x}).Excitation;
            INxFRFixed{x} = data.(IN{x}).Spikes;
            INxEXFree{x} = data_KCC2.(IN{x}).Excitation;
            INxFRFree{x} = data_KCC2.(IN{x}).Spikes;
            INxEX = {INxEXFixed, INxEXFree};
            if(not(and(data_KCC2.(IN{x}).Excitation(1)==0,data_KCC2.(IN{x}).Spikes(1)>5)))
                % Filter out excitatory inhibition
                INxFRFixed_filtered{x_index} = data.(IN{x}).Spikes;
                INxFRFree_filtered{x_index} = data_KCC2.(IN{x}).Spikes;
                x_index=x_index+1;
            end
            INxEGABAFixed{x} = data.(IN{x}).ldend_EGABA_;
            INxEGABAFree{x} = data_KCC2.(IN{x}).ldend_EGABA_;
        end
        [Nindex, fig,subfigs] = NormalizationIndex(INxEX,IN,INxFRFixed,INxFRFree,...
                            'titlename',options.titlename,...
                            'units',options.units,...
                            'plottype','semilogx',...
                            'plot',1,...
                            'format',options.format,...
                            'special',options.special,...
                            'color',options.color); 
        for x = 1:length(INxEGABAFree)
            % find difference in EGABA 
            % base EGABA taken from lowest excitation and lowest inhibition
            % data point (x) EGABA's value should be almost the same if retrieved from
            % (1) or (end) index - different EX values.
            DEGABA(x) = INxEGABAFree{x}(1) - INxEGABAFree{1}(1);
        end
        % Figure with non-proportional x-axis of inhibition
%         NvEGABAfig = figure();
%         plotyy(1:length(IN),Nindex,1:length(IN),DEGABA);
%         axis_right = NvEGABAfig.Children(1);
%         if (max(DEGABA)>0)
%             axis_right.YLim = [0, max(DEGABA)];
%         end
%         axis_left = NvEGABAfig.Children(2);
%         ax = gca;
%         ax.XTick = 1:1:length(IN);
%         ax.XTickLabels = multiReplaceAll(IN,{'c','d'},{'',''});
         result = Nindex;
%         if(and(toSave,toPlot))
%            saveFigs('folderName',{'NvEGABAfig'},(NvEGABAfig));
%         end
end

function [result,DEGABA_fig_number,NvEGABA_fig_number] = NormalisationIndexAnalysis(foldername, keepFunctions, varargin)
    global toSave toPlot
    warning('off', 'MATLAB:legend:IgnoringExtraEntries');
    defaultOptions = struct('static_cli',0,...
                            'heading','',...
                            'excludeSplit',{''},...
                            'excludeDegree',{''},...
                            'NindexInput',NaN,...
                            'plot',0);
    fixed_foldername = (multiReplace(foldername,{'_'},{''}));
    options = processOptions(defaultOptions,varargin);
    disp(fixed_foldername);
    colors = getColorPalette();
    brightens = [-0.6,-0.3,0.3,0.6];
    blue = palette('bluelight');
    red = palette('red');
    marker = 'o+sx*d^v><ph';
    line_styles = {'-','--',':','-.'};
    line_styles_len = length(line_styles);
    color_index = 1;
    marker_index = 1;
    plot_fittings=1;
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    % keep functions must be noKCC2 then KCC2
    for k = 1:2:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,...
            'excludeSplit',options.excludeSplit,...
            'excludeDegree',options.excludeDegree,...
            'trimBeginSplit',trimBeginSplit,...
            'trimEndSplit',trimEndSplit,...
            'trimBeginDegree',trimBeginDegree,...
            'trimEndDegree',trimEndDegree);
        if(options.static_cli>0)
           % static cli
            output_kcc2 = output;
            base_split_field = strcat('s',num2str(options.static_cli));
            try
                output.(base_split_field);
            catch
                error('%s not present.\n Found: %s',base_split_field,strjoin(fields(output)));
            end
        else
            keep = keeps{k+1};
            keepname = keepnames{k+1};
            output_kcc2 = getCompare(foldername,keep,splitReg,degree,...
                'excludeSplit',options.excludeSplit,...
                'excludeDegree',options.excludeDegree,...
                'trimBeginSplit',trimBeginSplit,...
                'trimEndSplit',trimEndSplit,...
                'trimBeginDegree',trimBeginDegree,...
                'trimEndDegree',trimEndDegree);
        end
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        [split_fields_kcc2, split_fields_num_kcc2] = reorderByNum(fields(output_kcc2));
        if (length(fields(output))==1 && length(fields(output_kcc2))>1)
            for f = 1:length(split_fields_kcc2)
                output.(split_fields_kcc2{f}) = output.(split_fields{1});
            end
        end
        
        [split_fields, split_fields_num] = reorderByNum(fields(output_kcc2));
        
        % group colors by limiting length of linestyles
        if strcmp(split_fields_num{1},split_fields_num{2})==1
            for i=2:length(split_fields_num)
                if strcmp(split_fields_num{1},split_fields_num{i}) ~= 1
                    line_styles_len = i-1;
                    break
                end
            end
        end
        maxDEGABA = 0.0000000001;
        DEGABAfig = figure();
        NvEGABAfig = figure();
        % get figure index to recall
        DEGABA_fig_number((k+1)/2) = DEGABAfig.Number;
        NvEGABA_fig_number((k+1)/2) = NvEGABAfig.Number;
        for split_index = 1:length(split_fields)
            if(options.static_cli>0)
                % static cli
                split_result = output.(base_split_field);                        
            else
                split_result = output.(split_fields{split_index});
            end
            split_result_kcc2 = output_kcc2.(split_fields{split_index});

            cutEnd=0;    
            [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
            for degree_index = 1:length(degree_fields_sub)
                %columns: same degree, different split (:,df)
                %rows: same split, different degree (rf,:)
                excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
                spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
                somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
                bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA;
                ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA_;

                excitation_kcc2{split_index,degree_index} = split_result_kcc2.(degree_fields_sub{degree_index}).Excitation; 
                spikes_kcc2{split_index,degree_index} = split_result_kcc2.(degree_fields_sub{degree_index}).Spikes;
                somaEGABA_kcc2{split_index,degree_index} = split_result_kcc2.(degree_fields_sub{degree_index}).Soma_EGABA;
                bdendEGABA_kcc2{split_index,degree_index} = split_result_kcc2.(degree_fields_sub{degree_index}).bdend_EGABA;
                ldendEGABA_kcc2{split_index,degree_index} = split_result_kcc2.(degree_fields_sub{degree_index}).ldend_EGABA_;
            end
            degree_fields{split_index} = degree_fields_sub;
            degree_fields_num{split_index} = degree_fields_sub_num;
            distal = split_result;
            distal_KCC2 = split_result_kcc2;
            IN = degree_fields_sub;

            x_index=1;
            for x=1:length(IN)
                INxEXFixed{x} = distal.(IN{x}).Excitation;
                INxFRFixed{x} = distal.(IN{x}).Spikes;
                INxEXFree{x} = distal_KCC2.(IN{x}).Excitation;
                INxFRFree{x} = distal_KCC2.(IN{x}).Spikes;
                if(not(and(distal_KCC2.(IN{x}).Excitation(1)==0,distal_KCC2.(IN{x}).Spikes(1)>5)))
                    % Filter out excitatory inhibition
                    INxFRFixed_filtered{x_index} = distal.(IN{x}).Spikes;
                    INxFRFree_filtered{x_index} = distal_KCC2.(IN{x}).Spikes;
                    x_index=x_index+1;
                end
                INxEGABAFixed{x} = distal.(IN{x}).ldend_EGABA_;
                INxEGABAFree{x} = distal_KCC2.(IN{x}).ldend_EGABA_;
            end
            EX = {INxEXFixed, INxEXFree};
            [Nindex, Nfig, Nsubfigs] = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree,...
               'titlename',split_fields_num{split_index},... 
               'units','times baseline conductance',...
               'plottype','semilogx',...
               'plot',options.plot); 
%            combinedfig{split_index} = combineFigures(figs,'titlename','',...
%                                         'plottype','logx',...
%                                         'matchcolors',0,...
%                                         'uniquespecifiers',1,...
%                                         'lines',1,...
%                                         'legend','full');
            
            if isstruct(options.NindexInput)
                try
                    replaceNindices = fields(options.NindexInput.(split_fields{split_index}));
                    for replaceNindex = replaceNindices
                       index = find(ismember(IN, replaceNindex)); 
                       Nindex(index) = options.NindexInput.(split_fields{split_index}).(replaceNindex{1});
                    end
                catch Error
                    %ignore error
                end
            end
            for x = 1:length(INxEGABAFree)
                % find difference in EGABA 
                % base EGABA taken from lowest excitation and lowest inhibition
                % data point (x) EGABA's value should be almost the same if retrieved from
                % (1) or (end) index - different EX values.
                if(options.static_cli>0)
                    %                               base cli, first inhibitory strength, first
                    %                               value of EGABA (they should all be the same
                    %                               for all strengths)
                    DEGABA(x) = INxEGABAFree{x}(1) - output.(base_split_field).(IN{1}).ldend_EGABA_(1);
                else
                    DEGABA(x) = INxEGABAFree{x}(1) - INxEGABAFree{1}(1);
                end
                DEGABA_cell(split_index,x) = DEGABA(x);
            end
            
            f_save(split_index) = figure();
            hold on;
            box off;
            
            x = toNumber(IN);
            IN_label = num2str(x);
            y = Nindex;
            % get rid of excitatory inhibition, NaNs, and Infs
            x = x(Nindex>=0);
            DEGABA = DEGABA(Nindex>=0);
            y = y(Nindex>=0);
            
            %% plot SMOOTH
%             y(Nindex<0) = [];
%             smoothed_y = smooth(x,y);
%             fitted_y = lineOfBestFit(x,y,2);
% 
%             plots_for_legend(1) = plot(x,y,'x','Color',colors(color_index,:));
%             if(plot_fittings)
%                 plots_for_legend(2) = plot(x,fitted_y,'-','Color',colors(color_index,:));
%             end
%             title(strcat(options.heading,'=',split_fields_num{split_index}));
%             % xlim([5,100]);
%             xlabel('Inhibition');
%             ylim([0,1]);
%             ylabel('Chloride Index');
%             % include fitted line in plots_for_legend to have values in legend and have
%             % it plotted in the combined figure
%             if(plot_fittings)
%                 legend(plots_for_legend,{split_fields_num{split_index},'(fitted)'});
%             else
%                 legend(plots_for_legend,split_fields_num(split_index));
%             end
            
            %% plot normalisation and delta EGABA
            yyFig = figure(DEGABA_fig_number((k+1)/2));
            edge_color = 'none';
            hold on;
            box off;
%             [AX,H1,H2] = plotyy(x,y,x,DEGABA);
%             H1.LineWidth = 2;
%             H2.LineWidth = 2;
%             H1.LineStyle = line_styles{color_index};
%             H2.LineStyle = line_styles{color_index};
%             H1.Color = blue;
%             H2.Color= red;
%             set(AX,{'ycolor'},{blue;red})  % Left color blue, right color red...
%             
            yyaxis right
            if color_index > line_styles_len
                marker_index = floor(color_index/line_styles_len);
            end
            color_r = hex2rgb('#D85218');
            color_b = hex2rgb('#0071BC');
            if strcmp(split_fields_num{1},split_fields_num{2})==1 
               color_r = brighten(color_r, brightens(marker_index));
               color_b = brighten(color_b, brightens(marker_index));
            end
            line_index = mod(color_index,line_styles_len)+1;
            plot(x,DEGABA, line_styles{line_index},...
                'Color',color_r,...
                'LineWidth',2);
            
            yyaxis left
            if length(x) ~= length(y)
                disp('error');
            end
            plot(x,y, line_styles{line_index},...
                'Color',color_b,...
                'LineWidth',2);
%             marker(color_index),
%             'MarkerSize',10,...
%                 'MarkerEdgeColor',edge_color,...
%                 'MarkerFaceColor',edge_color)
            
            if(max(DEGABA)>maxDEGABA)
                maxDEGABA = max(DEGABA); 
            end
            
            NvEGABAfig = figure(NvEGABA_fig_number((k+1)/2));
            hold on;
            box off;
            mdl = fitlm(DEGABA,y, 'Intercept', false);
            lines = plot(mdl);
            % return arguments for plot(mdl) are 4 lines
            % 1: Data
            % 2: Fit
            % 3: confidence bound bottom
            % 4: confidence bound top
            
            lines(1).Color = 'none';
            lines(2).Color = colors(color_index,:);
            lines(3).Color = 'none';
            lines(4).Color = 'none';
            
            %plot each point with a different marker to indicate different inhibitiory
            %strengths
            mi=1;
            data_markers = lines(1);
            for d=1:length(DEGABA)
                plot(DEGABA(d),Nindex(d),'Color', colors(color_index,:),...
                    'Marker',marker(mi),...
                    'MarkerFaceColor', colors(color_index,:));
                mi=mi+1;
                if(mi>length(marker))
                    mi=1;
                end
            end
            NvEGABAplots(split_index) = lines(2);
            
            color_index=color_index+1;
            if(color_index>length(colors))
                color_index=1;
            end
            
            result{(k+1)/2,split_index} = Nindex;
            N_mat(split_index,:) = Nindex;
        end
        
        fig((k+1)/2) = combineFigures(f_save,'titlename',options.heading,...
                                        'plottype','normal',...
                                        'matchcolors',0,...
                                        'uniquespecifiers',1,...
                                        'lines',1-plot_fittings,...
                                        'legend','full');
        if(plot_fittings)
            % we know that there are 2 lines per plot (the actual points and 
            % the fitted points). To make it look better, we shall change the 
            % fitted points to be a line ('lines' varargin was 0 above), and 
            % remove its marker
            lines = findobj(fig,'type','line');
            for l = 1:2:length(lines)
               lines(l+1).LineStyle = '-';
               lines(l+1).Marker = 'none';
               lines(l).LineStyle = ':';  % create a very faint line
            end
        end
        ylim([0,1]);
        xlabel('Inhibition');
        
        % add titles, labels, legends to DEGABA_fig
        yyFig = figure(DEGABA_fig_number((k+1)/2));
        title(options.heading);
        yyaxis right;
        set(gca,'linewidth',1)
        ylabel('\Delta EGABA (mV)');
        %ylim([0,50]);
        yyaxis left;
        %ylim([0,0.7]);
        ylabel('Chloride Index');
        xlabel('Relative inhibitory conductance');
        if str2double(split_fields_num(1))<0.01
            format = '.0e'; % round
        else
            format = '.2f';
        end
        if length(split_fields_num)>1 && strcmp(split_fields_num{1},split_fields_num{2})==1
           
            for s=1:length(split_fields)
                nums = split('c',split_fields{s}(2:end));
                for n=1:length(nums)
                    nums{n} = num2str(toNumber(nums{n}),strcat('%',format));
                end
                number{s} = strjoin(nums,',');
            end
            
        else
            number = strsplit(num2str(str2double(split_fields_num),strcat('%',format,'\n')),'\n');
        end
        legend(number);
        
        % add titles, labels, legends to NvEGABAfig
        NvEGABAfig = figure(NvEGABA_fig_number((k+1)/2));
        hold on;
        box off;
        set(gca,'linewidth',1)
        title(options.heading);
        xlabel('\Delta EGABA (mV)');
        ylabel('Chloride Index');
        ylim([0,1]);
        XL = xlim;
        legend(NvEGABAplots, number);
        %xlim([0,XL(2)]);
        %leg_main=legend(NvEGABAplots,split_fields_num,'location','northwest');
        %title(leg_main,options.heading)
        % add markings to indicate common inhibition for each DEGABA point
        % with connecting line: traces=plot(DEGABA_cell,N_mat,'x-');
        %traces=plot(DEGABA_cell,N_mat,':k','LineWidth',0.5);
        result{'d'} = DEGABA_cell;
        result{'n'} = N_mat;
        
        mdl = fitlm(DEGABA_cell(:),N_mat(:), 'Intercept', false);
        disp(mdl);
        text(5,0.9,strcat('R^2 = ',num2str(mdl.Rsquared.Ordinary)));
        
        % create new axis, and hide it
        ah=axes('position',get(gca,'position'));
        set(ah,'visible','off');
        %leg_inh = legend(ah,traces,IN_label','location','west');
        %title(leg_inh,'Inhibition');
   
        
        
        % save
        if(and(toSave,toPlot))
           %saveFigs('folderName',split_fields,f);
           f_save = figure(DEGABA_fig_number((k+1)/2));
           saveFigs('folderName',{strcat('DEGABA',options.heading,num2str((k+1)/2))},f_save);
           
           f_save = figure(NvEGABA_fig_number((k+1)/2));
           saveFigs('folderName',{strcat('NvEGABA',options.heading,num2str((k+1)/2))},f_save);
           ylim('auto');
           saveFigs('folderName',{strcat('NvEGABA',options.heading,num2str((k+1)/2),'auto')},f_save);
           
           f_save = fig((k+1)/2);
           saveFigs('folderName',{keepname},f_save);
        end
    end
    
    warning('on', 'MATLAB:legend:IgnoringExtraEntries');
end

function [result, fig] = Synapse_Numbers(foldername, keepFunctions, colorPalette, varargin)
    global folderName
    folderName = foldername;
    defaultOptions = struct('excludeInh',{''},...
                            'excludeExc',{''});
    options = processOptions(defaultOptions,varargin);
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();    
    inh_hz_index = 5;
    exc_hz_index = 5;
    inh_hz_index_str = strcat('d',num2str(inh_hz_index));
    
    line_spec={'--','-'};
    for k = 1:length(keeps)
        f = figure();
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,...
            'trimBeginSplit',trimBeginSplit,...
            'trimEndSplit',trimEndSplit,...
            'trimBeginDegree',trimBeginDegree,...
            'trimEndDegree',trimEndDegree);
        % transpose data so that strength is by synapse number instead of
        % frequency
        output = transposeSynapticNumberToStrength(output);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        slope=0;
        for split_index = 1:length(split_fields)
            if (strcmp(split_fields{split_index},strcat('e',num2str(exc_hz_index),'i',num2str(inh_hz_index))))
                top_split_index = split_index;
                split_result = output.(split_fields{split_index});
                [degree_fields, degree_fields_num] = reorderByNum(fields(split_result));
                
                for r=1:length(options.excludeInh)
                    % see if exclude in degree_fields_num
                    indexC = strfind(degree_fields_num, options.excludeInh{r});
                    % get indices
                    index = find(not(cellfun('isempty',indexC)));
                    if ~isempty(index)
                        % remove values
                        index=index(1);
                        degree_fields_num(index) = [];
                        degree_fields(index) = [];
                    end
                end
                
                colors = getColorPalette('Palette',colorPalette,...
                                         'Black',true, ...
                                         'N', length(fields(split_result)));
                c=0;
                for degree_index = 1:length(degree_fields)
                    degree_result = split_result.(degree_fields{degree_index});
                    excitation = degree_result.Excitation;
                    spikes = degree_result.Spikes;
                    [argvalue, argmax] = max(spikes);
                    if (argvalue <= 1)
                       argmax=length(spikes);
                    end
                    if (argmax < length(spikes))
                        if not(spikes(argmax) >= spikes(end)+5)
                            % max reached and sustained (within 5 hz)
                            argmax = length(spikes);
                        end
                    end
                    spikes = spikes(1:argmax);
                    excitation = excitation(1:argmax);
                    top_excitation{degree_index,k} = excitation;
                    top_spikes{degree_index,k} = spikes;
                end
            end
        end
        split_fields = multiReplaceAll(split_fields,{'e','i'},{'E:',' I:'}); 
    end
    
   % Calculate Normalization Index and plot combined figure
    IN = degree_fields;
    % which is fixed and which is free based on order of keepfunctions
    INxEXFixed = top_excitation(:,1);
    INxEXFree = top_excitation(:,2);
    INxFRFixed = top_spikes(:,1);
    INxFRFree = top_spikes(:,2);
    EX = {INxEXFixed, INxEXFree};
    IN2 = {INxFRFixed, INxFRFree};
    [a,b] = reorderByNum(degree_fields);
    [Nindex, fig] = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree,...
        'units','# synapses',...
        'format','%6.0f',...
        'color',colorPalette,...
        'synapses',1);
    result = Nindex;
end
function result = Balanced_input(foldername, inh_loc, varargin)
    global folderName toSave
    folderName = foldername;
    
    defaultOptions = struct('plotUnbalanced',0,'freq',20,...
        'excludeSplit',{''},'excludeDegree',{''});
    options = processOptions(defaultOptions,varargin);
    
    if strcmp(inh_loc, 'Proximal') || strcmp(inh_loc, 'proximal')
        keepFunctions = {@getKeepproximalKCC2,@getKeepproximalNoKCC2};
        colorPalette = 'Greens';
    elseif strcmp(inh_loc, 'Distal')
        keepFunctions = {@getKeepDistalKCC2,@getKeepDistalNoKCC2};
        colorPalette = 'Blues';
    else
        result = -1;
        return
    end
    
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    linear = 1:100;
    off_linf = figure();
    hold on;
    linestyles={'-','--'};

    % GET DATA and do some prelim plotting
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,...
            'excludeSplit',options.excludeSplit,...
            'excludeDegree',options.excludeDegree,...
            'trimBeginSplit',trimBeginSplit,...
            'trimEndSplit',trimEndSplit,...
            'trimBeginDegree',trimBeginDegree,...
            'trimEndDegree',trimEndDegree);
        topoutput.(funcnames{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           colors = cbrewer('seq', colorPalette, length(split_fields)+1);
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);

           if options.plotUnbalanced
               if k==1
                   unbalancedFigs{split_index} = figure();
                   hold on;
                   title(multiReplaceAll(split_fields{split_index},{'s','c'},{'E:',' I:'}));
               else
                   figure(unbalancedFigs{split_index}.Number);
               end
           end
           
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               degree_result = split_result.(degree_fields_sub{degree_index});
               exc = split_result.(degree_fields_sub{degree_index}).Excitation;
               [exc, indices] = sort(exc);
               inh =  str2num(degree_fields_sub_num{degree_index});
               spks = degree_result.Spikes;
               spks = spks(indices);
               % get balanced input value (exc == inh)
               [argvalue, argmin] = min(abs(inh - exc));
               
               excitation(split_index,degree_index) = exc(argmin); 
               spikes(split_index,degree_index) = spks(argmin);
               
               if options.plotUnbalanced
                   plot(exc,spks, ...
                       'LineStyle',linestyles{k},...
                       'LineWidth',2,...
                       'Color', [colors(split_index+1,:), 0.2])
               end
               
           end
           if options.plotUnbalanced
               p = plot(excitation(split_index,:), spikes(split_index,:),...
                   'LineStyle',linestyles{k},...
                   'LineWidth',2,...
                   'Color',colors(split_index+1,:));
               if k==1
                legend([p],['balanced']);
               end
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           
           off_linear_exc = excitation(split_index,:) - 0:5:5*length(excitation(1,:));
           lin = excitation(split_index,:);
           off_linear_spikes = spikes(split_index,:) - lin;
           figure(off_linf.Number);
           plot(excitation(split_index,:),off_linear_spikes,...
               'LineStyle',linestyles{k},...
               'LineWidth',2,...
               'Color',colors(split_index+1,:))
        end
        topoutput.(funcnames{k}).excitation = excitation;
        topoutput.(funcnames{k}).spikes = spikes;
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'',':'});
        for i=1:length(split_fields)
            colon_index = find(split_fields{i}==':');
            if(length(split_fields{i}(colon_index+1:end))==2)
                split_fields{i} = sprintf('%3s% 3s',...
                    split_fields{i}(1:colon_index),...
                    split_fields{i}(colon_index+1:end));
            end
        end
        %legend(split_fields);
        xlabel('Balanced Input (Hz)');
        ylabel('Output - Input (Hz)');
        
        result = 1;
    end
    title('Off-linear');
    off_linf=formatFigure(off_linf);
    
    
    figure(off_linf);
    plot(0:5:5*length(excitation(1,:))-1,linspace(0,0,length(spikes(1,:))),...
        'k-','LineWidth',2);
   
    new_off_linf = formatFigure(off_linf);
    new_off_linf.CurrentAxes.XLim = [0,60];
    
    if(toSave)
        saveFigs('folderName',{strcat('offlinear_',inh_loc)},off_linf);
    end
    
    l = findobj(new_off_linf,'type','line');
    leg = findobj(new_off_linf,'type','legend');
    colors = cbrewer('seq',colorPalette, length(split_fields)+1);
    
%     for split_index = 1:length(split_fields)
%         f2(split_index) = figure();
%         hold on;
%         for k=1:length(keeps)
%            plot(topoutput.(funcnames{k}).excitation(split_index,:),topoutput.(funcnames{k}).spikes(split_index,:))
%         end
%         plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
%         legend(keepnames);
%         xlabel('Input (Hz)');
%         ylabel('Output (Hz)');
%         title(split_fields{split_index});
%     end
%     saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(funcnames{1}).excitation(split_index,:),topoutput.(funcnames{1}).spikes(split_index,:),'-','Color',colors(split_index,:),'LineWidth',2);
        plot(topoutput.(funcnames{2}).excitation(split_index,:),topoutput.(funcnames{2}).spikes(split_index,:),'--','Color',colors(split_index,:),...
            'LineWidth',2);
    end
%     f3.CurrentAxes.XLim = [0,60];
%     f3.CurrentAxes.YLim = [0,30];

    % plot linear line
    plot(0:5:excitation(1,end),0:5:excitation(1,end),'k-',...
        'LineWidth',2);
    leg = legend(p,split_fields);
    title(leg,{'# Synapses';sprintf('\t(E:I)')})
    xlabel('Balanced Input (Hz)');
    ylabel('Output (Hz)');
    title(inh_loc);
    xlim([0, 60]);
    f3 = formatFigure(f3);
    if(toSave)
        saveFigs('folderName',{strcat('balancedio_',inh_loc)},f3);
    end
    
    % plot SATURATION
    f4 = figure();
    hold on;
    freq = options.freq;
    for split_index = 1:length(split_fields)
        exc = topoutput.(funcnames{2}).excitation(split_index,:);
        spikes = topoutput.(funcnames{2}).spikes(split_index,:);
        exc_KCC2 = topoutput.(funcnames{1}).excitation(split_index,:);
        spikes_KCC2 = topoutput.(funcnames{1}).spikes(split_index,:);
        indexof = find(exc==freq);
        inh_synapses(split_index) = str2double(split_fields{split_index}(find(split_fields{split_index}==':')+1:end));
        inh_syn_spikes(split_index) = spikes(indexof);
        inh_syn_spikes_KCC2(split_index) = spikes_KCC2(indexof);
    end
    %plot(inh_synapses,inh_syn_spikes,'--',inh_synapses,inh_syn_spikes_KCC2,'-');
    f = fit(inh_synapses',inh_syn_spikes','exp2');
    f_KCC2 = fit(inh_synapses',inh_syn_spikes_KCC2','exp2');
    p = plot(f,'k--',inh_synapses,inh_syn_spikes,'.');
    p_KCC2 = plot(f_KCC2,'k-',inh_synapses,inh_syn_spikes_KCC2,'.');
    
    % change color of lines and fill in circles
%     p(2).Color=[0.5 0.5 0.5];
%     p_KCC2(2).Color=[0.5 0.5 0.5];
    for i=1:length(inh_synapses)
       plot( inh_synapses(i),inh_syn_spikes(i),'.',...
           'color',colors(i+1,:),'MarkerFaceColor',colors(i+1,:),'MarkerSize',35.0);
       plot( inh_synapses(i),inh_syn_spikes_KCC2(i),'.',...
           'color',colors(i+1,:),'MarkerFaceColor',colors(i+1,:),'MarkerSize',35.0)
    end
    p(1).MarkerFaceColor='none';
    
    p_KCC2(1).MarkerFaceColor='k';
    
    leg = legend([p_KCC2(2),p(2)],{'Dynamic Chloride','Static Chloride'});
    xlabel('Inhibitory synapses (#)');
    ylabel(strcat('Output given ', num2str(freq), ' Hz input (Hz)'));
    title(strcat('Input: ', num2str(freq), ' Hz'));
    formatFigure(f4);
    if(toSave)
        saveFigs('folderName',{strcat('saturation_',inh_loc,'_',options.freq)},f4);
    end
    % plot OUTPUT RANGE
    f5 = figure();
    hold on;
    for split_index = 1:length(split_fields)
        exc = topoutput.(funcnames{2}).excitation(split_index,:);
        spikes = topoutput.(funcnames{2}).spikes(split_index,:);
        exc_KCC2 = topoutput.(funcnames{1}).excitation(split_index,:);
        spikes_KCC2 = topoutput.(funcnames{1}).spikes(split_index,:);
        max_plot(split_index,1) = max(spikes);
        max_plot(split_index,2) = max(spikes_KCC2)-max(spikes);
        
    end
    b = bar(inh_synapses,max_plot,'stacked');
    b(1).LineStyle = '--';
    b(2).FaceColor = [0,0,0];
    b(1).FaceColor = [1,1,1];
    b(1).LineWidth = 2;
    b(2).LineWidth = 2;
%     f5.CurrentAxes.XLim = [0,160];
    %xticklabels(split_fields)
    %set(gca,'XTickLabel', split_fields)
    %leg = legend([p(2),p_KCC2(2)],{'No','Yes'});
    title(leg,'Chloride dynamics');
    xlabel('Inhibitory synapses (#)');
    ylabel('Maximum Output (Hz)');
    title(strcat('output range'));
    formatFigure(f5);
    if(toSave)
        saveFigs('folderName',{strcat('maxoutput_',inh_loc)},f5);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REPEATABLE DATA RUNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = proximal_balanced_synapses_figure()
    % run pattern of inhibition for matched input, but using
    %  optimal synaptic numbers
    whos
    foldername = 'Hz/Data_balanced_synaptic_input_proximal';
    result = Balanced_input(foldername, 'Proximal', 'plotUnbalanced', 0,...
        'freq',50,...
        'excludeSplit',{'370,120','430,150','530,210','600,240','720,300','660,270'});
end
function result = distal_balanced_synapses_figure()
    % run pattern of inhibition for matched input, but using
    %  optimal synaptic numbers
    whos
    foldername = 'Hz/Data_balanced_synaptic_input_distal';
    result = Balanced_input(foldername, 'Distal', 'plotUnbalanced', 0,...
        'excludeSplit',{'220,90','240,200','250,260','260,400','280,600','290,700'},...
        'excludeDegree',{''});
end

function [result,fig] = Data_persistent_somatic()   
    whos
    foldername = 'Hz/Data_2019_11_27_16h18_persistent';
    keepFunctions = {@getKeepSomaticNoKCC2};   
     for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();    
    fig = figure();
    keep = keeps{1};
    keepname = keepnames{1};
    output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
    % transpose data so that strength is by synapse number instead of
    % frequency
    %output = transposeSynapticNumberToStrength(output);
    [split_fields, split_fields_num] = reorderByNum(fields(output));

    slope=0;
    for split_index = 1:length(split_fields)
       
        top_split_index = split_index;
        split_result = output.(split_fields{split_index});
        [degree_fields, degree_fields_num] = reorderByNum(fields(split_result));

        colors = getColorPalette('Palette','Reds',...
                                 'Black',true, ...
                                 'N', length(fields(split_result)));
        c=0;
        for degree_index = 1:length(fields(split_result))
            degree_result = split_result.(degree_fields{degree_index});
            excitation = degree_result.Excitation;
            spikes = degree_result.Spikes_at_Soma;
            [argvalue, argmax] = max(spikes);
            if (argvalue <= 1)
               argmax=length(spikes);
            end
            spikes = spikes(1:argmax);
            excitation = excitation(1:argmax);

            semilogx(excitation,spikes,'-','Marker','o');
            hold on;
            box off;
        end
        
    end
    split_fields = multiReplaceAll(split_fields,{'e','i'},{'E:',' I:'}); 
    
   % Calculate Normalization Index and plot combined figure
    IN = degree_fields;
    legend(IN);
end
function [result,fig] = Data_Increasing_synapse_numbers_distal()
    % 5 Hz input.
    % increase numbers of synapses to determine gain shift
    % proximal
    % with + without KCC2
    % synapses are evenly distributed
    % weights are (E:I) = (2:7)
    
    whos
    foldername = 'Hz/Data_Increasing_synapse_numbers';
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [result, fig] = Synapse_Numbers(foldername, keepFunctions, 'Blues',...
                        'excludeInh',{'10','20','40','158'});
    fig = formatFigure(fig);
    
end
function [result,fig] = Data_Increasing_synapse_numbers_proximal()
    % 5 Hz input.
    % increase numbers of synapses to determine gain shift
    % proximal
    % with + without KCC2
    % synapses are evenly distributed
    whos
    foldername = 'Hz/Data_Increasing_synapse_numbers';
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    [result, fig] = Synapse_Numbers(foldername, keepFunctions, 'Greens',...
                'excludeInh',{'10','20','158','630','998'});

end

function result = n_index_figure_traces()
    % create traces from txt file
    global folderName
    whos
    foldername = 'Traces';
    folderName = foldername;
    keepFunctions = {@getKeepDistal,@getKeepProximal};
    [split,splitName,trimBeginSplit,trimEndSplit]=getSplitOnKCC2();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        disp(output);
        f1 = figure();
        plot(output.sno.d.time,output.sno.d.voltage,'k');
        title(strcat(keepname,'noKCC2'));
        xlim([0,1000]);
        formatFigure(f1);
        f2 = figure();
        plot(output.sKCC2.d.time,output.sKCC2.d.voltage,'k');
        title(strcat(keepname,'KCC2'));
        xlim([0,1000]);
        formatFigure(f2);
    end
    result =1;
end
function result = n_index_figure()
    global toSave
    % create figures based on the protocols below
    % to be combined in vector editing software
    excludeINdistal = {'2.5','7.5','12.5','15','17.5','25','100','200','500'}; 
    highlight={'20'};
    debug_plot = 0; % 0 for none, 1 for NIndex figures, 2 for individual IN comparison
    result=1;
    % DISTAL N-index
   %[Nindex,nIndex_fig,subfigs] = Data_persistent_distal('excludeDegree',excludeINdistal,'special',highlight);
    % PROXIMAL N-index
    excludeINdistal_prox = {'0.005','0.01','0.05','5.0','17.5','100','200','500'}; 
    %Data_persistent_proximal('excludeDegree',excludeINdistal_prox, 'special',{'0.5'});
    %specialNindex = struct('s1',struct('d5',Nindex(2:end)));
    specialNindex= NaN;
    NvDEGABA_nums = [];  % note these are figure numbers
    results_degaba = [];
    results_nidx = [];
    functions = {@Data_persistent_pkcc2, @Data_persistent_diam,...
        @Data_persistent_pas,@Data_persistent_pkcc2_homo,...
        @Data_persistent_duration, @Data_persistent_dynamic_K};
    others = {@Data_persistent_pas,@Data_persistent_pkcc2_homo,...
        @Data_persistent_duration,@Data_persistent_dynamic_K};
    names = {'KCC2 strength', 'Distal diameter', 'IR', 'KCC2 homogeneity', ...
        'Duration', 'K^+'};
    %Data_2018_02_02_23h48_persistent_pcl Data_2018_02_03_02h31_persistent_cli
    
    num_lines = [];
    f_combined = figure();
    hold on;
    colors = getColorPalette();
    specifier = 'osd+x*.^v><ph';
    
    for f = 1:length(functions)
        [result,DEGABA_num,NvDEGABA_num] = functions{f}({'100','200','500'}, specialNindex, debug_plot);
        figure(NvDEGABA_num);
        NvDEGABA_nums = [NvDEGABA_nums;NvDEGABA_num];
        dgaba = result{'d'};
        nidx= result{'n'};
        flat_dgaba = reshape(dgaba,1,numel(dgaba));
        flat_nidx = reshape(nidx,1, numel(nidx));
        results_degaba = [results_degaba flat_dgaba];
        results_nidx = [results_nidx flat_nidx];
        lines = findobj(NvDEGABA_num,'type','line');
        num_lines = [num_lines; length(lines)];
        figure(f_combined.Number);
        plot(flat_dgaba, flat_nidx, specifier(f), ...
            'Color', colors(f,:), ...
            'MarkerFaceColor', colors(f,:));
    end
    
    % perform stats
    obs = numel(results_nidx);
    [R,P] = corrcoef(results_degaba,results_nidx);
    stats = strcat('r=',num2str(R(2)),', p=',num2str(P(2)),', Pearson correlation, N=', num2str(obs));
    disp(stats);
    
    mdl = fitlm(results_degaba(:),results_nidx(:));
    disp(mdl);
    
    % plot combined NIndex vs DEGABA
    
    legend(names{1:length(functions)});
    ax = findobj(NvDEGABA_num,'type','axes');
    xlabel(ax(2).XLabel.String);
    ylabel(ax(2).YLabel.String);
    text(0,0.4,stats)
    
    if(toSave)
        saveFigs('folderName',{'NvDEGABA_combined'},(f_combined));
    end
    
end
function [result,fig,subfigs] = Data_persistent_distal(varargin)
    % use persistent synapses as a proxy for number of synapses
    % better values for excitation and inhibition strength
    global folderName toPlot toSave
    whos
    defaultOptions = struct('excludeSplit',{''},...
                            'excludeDegree',{''},...
                            'special',{''});
    options = processOptions(defaultOptions,varargin);
    foldername = 'Hz/Data_persistent';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    %% KCC2 Comparison
    prox = 0;
    dist = 1;
    inter = 0;
    toplot = 0;
    save = 0;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save,'semilogx',...
        'excludeSplit',options.excludeSplit,'excludeDegree',options.excludeDegree);
    
    if(dist)
        IN = fields(distal);
        EX = distal.(IN{1}).Excitation;
        [result,fig,subfigs] = NormalisationIndexAnalysisBase(distal,distal_KCC2,EX,IN,...
        'titlename','distal',...
        'format','%.1f',...
        'special',options.special,...
        'units','relative conductance');
    end
    
    
    
end
function [result,fig,subfigs] = Data_persistent_proximal(varargin)
    % use persistent synapses as a proxy for number of synapses
    % better values for excitation and inhibition strength
    global folderName toPlot toSave
    whos
    defaultOptions = struct('excludeSplit',{''},...
                            'excludeDegree',{''},...
                            'special',{''});
    options = processOptions(defaultOptions,varargin);
    foldername = 'Hz/Data_persistent';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    %% KCC2 Comparison
    prox = 1;
    dist = 0;
    inter = 0;
    toplot = 0;
    save = 0;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save,'semilogx',...
        'excludeSplit',options.excludeSplit,'excludeDegree',options.excludeDegree);
    
    if(prox)
        IN = fields(proximal);
        EX = proximal.(IN{1}).Excitation;
        [result,fig,subfigs] = NormalisationIndexAnalysisBase(proximal, proximal_KCC2,EX,IN,...
        'titlename','proximal',...
        'format','%.2f',...
        'color','Greens',...
        'special',options.special,...
        'units','relative conductance');
    end
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_pkcc2_homo(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_pkcc2_homo';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{},...
        'heading','KCC2 Homogeneity',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    disp(NvEGABA_num);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_duration(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_duration';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'100.0'},...
        'heading','duration (ms)',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_dynamic_K(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_K';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{},...
        'heading','dynamic K^+',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_pas(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_pas';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'[1.5,255.51,4.43,-70.75]',...
                        '[2,197.47,4.59,-70.51]',...
                        '[3,136.75,4.89,-70.10]',...
                        '[5,85.55,5.39,-69.46]',...
                        '[8,55.32,6.00,-68.76]',...
                        '[10,44.93,6.33,-68.40]'},...
        'heading','passive',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_pas_diam1(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_pas_diammix';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'[1.5,255.51,4.43,-70.75]',...
                        '[2,197.47,4.59,-70.51]',...
                        '[3,136.75,4.89,-70.10]',...
                        '[5,85.55,5.39,-69.46]',...
                        '[8,55.32,6.00,-68.76]',...
                        '[10,44.93,6.33,-68.40]'},...
        'heading','passive',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_pcl(excludeINdistal,specialNindex,debug_plot)
    %% Internal chloride ion concentration comparison
    % manually set for each run
    global folderName
    whos
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_pcl';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};                
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{''},...
        'heading','[Cl^-] permeability',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_cli(excludeINdistal,specialNindex,debug_plot)
    %% Internal chloride ion concentration comparison
    % manually set for each run
    global folderName
    whos
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_cli';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};                
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{},...
        'static_cli',5,...
        'heading','static [Cl^-]',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_pkcc2(excludeINdistal,specialNindex,debug_plot)
    %% pkcc2 comparison
    global folderName
    whos
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_pkcc2';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};                
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{''},...
        'heading','KCC2 strength',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num,NvEGABA_num] = Data_persistent_diam(excludeINdistal,specialNindex,debug_plot)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
        debug_plot=0;
    end
    foldername = 'Hz/Data_persistent_diam';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num,NvEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'0.25','0.4'},...
        'heading','distal diameter (um)',...
        'NindexInput',specialNindex,...
        'plot',debug_plot);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end

%% OTHER
function result = Data_2017_03_11_16h48_persistent()
    % use persistent synapses as a proxy for number of synapses
    % better values for excitation and inhibition strength
    global folderName toPlot toSave
    whos
    foldername = 'Hz/Data_2017_03_11_16h48_persistent';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    %% KCC2 Comparison
    prox = 0;
    dist = 1;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save,'semilogx');
    if(dist)
        IN = fields(distal);
        EX = distal.(IN{1}).Excitation;
        result = NormalisationIndexAnalysisBase(distal,distal_KCC2,EX,IN,'titlename','distal');
    end

end
function [result,DEGABA_num] = Data_2016_11_08_19h29_persistent_cli(excludeINdistal)
    %% Internal chloride ion concentration comparison
    % manually set for each run
    global folderName
    whos
    foldername = 'Hz/Data_2016_11_08_19h29_persistent_cli';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};                
    [result,DEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'2.5','20','25','30'},...
        'static_cli',5,...
        'heading','static [Cl^-]');
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num] = Data_2016_10_28_19h56_persistent_pkcc2(excludeINdistal)
    %% pkcc2 comparison
    global folderName
    whos
    % Data_2016_10_28_19h56_persistent_pkcc2
    foldername = 'Hz/Data_2016_10_28_19h56_persistent_pkcc2';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};                
    [result,DEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{''},...
        'heading','KCC2 strength');
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function [result,DEGABA_num] = Data_2016_10_25_19h01_persistent_diam(excludeINdistal,specialNindex)
    %% diam comparison
    global folderName
    whos 
    if(nargin==1)
        specialNindex = NaN;
    end
    foldername = 'Hz/Data_2016_10_25_19h01_persistent_diam';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    [result,DEGABA_num] = NormalisationIndexAnalysis(foldername,keepFunctions,...
        'excludeDegree',excludeINdistal,...
        'excludeSplit',{'0.25','0.4'},...
        'heading','distal diameter (um)',...
        'NindexInput',specialNindex);
    % no longer useful for analytics, but good for graphing
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
end
function result = Data_2016_10_24_21h21_persistent()
    % use persistent synapses as a proxy for number of synapses
    % better values for excitation and inhibition strength
    global folderName toPlot toSave
    whos
    foldername = 'Hz/Data_2016_10_24_21h21_persistent';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    %% KCC2 Comparison
    prox = 1;
    dist = 1;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save,'semilogx');
    % Get excitation and spikes (at axon)
    if(prox)
        IN = fields(proximal);
        EX = proximal.(IN{1}).Excitation;     %Excitation values are the same for all proximal runs
        result = NormalisationIndexAnalysisBase(proximal,proximal_KCC2,EX,IN,'titlename','proximal');
    end
    if(dist)
        IN = fields(distal);
        EX = distal.(IN{1}).Excitation;           %Excitation values are the same for all distal runs
        result = NormalisationIndexAnalysisBase(distal,distal_KCC2,EX,IN,'titlename','distal');
    end
%     for x=1:length(proximal_fields)
%         INx1{x} = proximal.(proximal_fields{x}).Spikes;
%         INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
%     end
%     
%     % Generate Nindex
%     [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
%     
%     % Check 'desired firing rate' curve
%     inhibitonX = 1;
%     [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
%     [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
   
end
function result = Data_2016_10_20_14h10_persistent()
    % use persistent synapses as a proxy for number of synapses
    global folderName toPlot toSave
    whos
    foldername = 'Hz/Data_2016_10_20_14h10_persistent';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    %% KCC2 Comparison
    prox = 0;
    dist = 1;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save,'semilogx');
    %cd(old);
    % Get excitation and spikes (at axon)
%     proximal_fields = fields(proximal);
%     EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        IN_distal = fields(distal);
        EX_distal = distal.(IN_distal{1}).Excitation;           %Excitation values are the same for all distal runs
    end
%     for x=1:length(proximal_fields)
%         INx1{x} = proximal.(proximal_fields{x}).Spikes;
%         INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
%     end
%     
%     % Generate Nindex
%     [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
%     
%     % Check 'desired firing rate' curve
%     inhibitonX = 1;
%     [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
%     [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    
    if(dist) 
        x_index=1;
        for x=1:length(IN_distal)
            INxFRFixed{x} = distal.(IN_distal{x}).Spikes;
            INxFRFree{x} = distal_KCC2.(IN_distal{x}).Spikes;
            if(not(and(distal_KCC2.(IN_distal{x}).Excitation(1)==0,distal_KCC2.(IN_distal{x}).Spikes(1)>5)))
                % Filter out excitatory inhibition
                INxFRFixed_filtered{x_index} = distal.(IN_distal{x}).Spikes;
                INxFRFree_filtered{x_index} = distal_KCC2.(IN_distal{x}).Spikes;
                x_index=x_index+1;
            end
            INxEGABAFixed{x} = distal.(IN_distal{x}).ldend_EGABA_;
            INxEGABAFree{x} = distal_KCC2.(IN_distal{x}).ldend_EGABA_;
        end
        [Nindex_distal, fig] = NormalizationIndex(EX_distal,IN_distal,INxFRFixed,INxFRFree); 
        for x = 1:length(INxEGABAFree)
            % find difference in EGABA 
            % base EGABA taken from lowest excitation and lowest inhibition
            % data point (x) EGABA's value should be almost the same if retrieved from
            % (1) or (end) index - different EX values.
            DEGABA(x) = INxEGABAFree{x}(1) - INxEGABAFree{1}(1);
        end
        f = figure();
        plotyy(1:length(IN_distal),Nindex_distal,1:length(IN_distal),DEGABA);
        axis_right = f.Children(1);
        axis_right.YLim = [0, max(DEGABA)];
        axis_left = f.Children(2);
        ax = gca;
        ax.XTick = 1:1:length(IN_distal);
        ax.XTickLabels = multiReplaceAll(IN_distal,{'c','d'},{'',''});
        result = Nindex_distal;
    end 
    
    
end
function result = Data_2016_10_12_14h31_protocol_increasing_weights()
    % 5 Hz input.
    % increase numbers of synapses (through weight value) to determine gain shift
    % proximal + distal
    % with + without KCC2
    % 10 initial synapses each and weight adjusted
    global folderName
    whos
    foldername = 'Hz/Data_2016_10_12_14h31_protocol_increasing_weights';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2,@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [split,splitName,trimBeginSplit,trimEndSplit]=getSplitOnWeights();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    inh_hz_index = 5;
    exc_hz_index = 5;
    inh_hz_index_str = strcat('d',num2str(inh_hz_index));
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        % transpose data so that strength is by synapse number instead of
        % frequency
        output = transposeSynapticNumberToStrength(output);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        f(k) = figure();
        
        for split_index = 1:length(split_fields)
            if (strcmp(split_fields{split_index},strcat('e',num2str(exc_hz_index),'i',num2str(inh_hz_index))))
                top_split_index = split_index;
                split_result = output.(split_fields{split_index});
                [degree_fields, degree_fields_num] = reorderByNum(fields(split_result));
                for degree_index = 1:length(fields(split_result))
                    degree_result = split_result.(degree_fields{degree_index});
                    excitation = degree_result.Excitation;
                    spikes = degree_result.Spikes;
                    top_excitation{degree_index,k} = excitation;
                    top_spikes{degree_index,k} = spikes;
                    semilogx(excitation,spikes);
                    hold on;
                end
                leg = legend(degree_fields_num);
            end
        end
        split_fields = multiReplaceAll(split_fields,{'e','i'},{'E:',' I:'}); 
        title(keepname)
        title(leg,['Inhibition' char(10) '(# channels per synapse)']);
        xlabel('Excitation (# channels per synapse)');
        ylabel('Output (Hz)');
    end
%     
%     % Calculate Normalization Index
%     EX = top_excitation{top_split_index,1};
%     IN = degree_fields;
%     % which is fixed and which is free based on order of keepfunctions
%     INxFRFixed = top_spikes(:,1);
%     INxFRFree = top_spikes(:,2);
%     Nindex = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree);
%     result = Nindex;
    proximalFig = combineFigures(f(1:2),...
                            'titlename','proximal',...
                            'plottype','logx',...
                            'matchcolors',1,...
                            'uniquespecifiers',0,...
                            'lines',1,...
                            'legend','once');
    distalFig = combineFigures(f(3:4),...
                            'titlename','distal',...
                            'plottype','logx',...
                            'matchcolors',1,...
                            'uniquespecifiers',0,...
                            'lines',1,...
                            'legend','once');
    result =1;
end
function result = Data_2016_10_07_12h38_Increasing_synapse_numbers()
    % 5 Hz input.
    % increase numbers of synapses to determine gain shift
    % distal
    % with + without KCC2
    % synapses are evenly distributed
    % weights are (E:I) = (2:7)
    global folderName toPlot
    whos
    foldername = 'Hz/Data_2016_10_07_12h38_Increasing_synapse_numbers';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    inh_hz_index = 5;
    exc_hz_index = 5;
    inh_hz_index_str = strcat('d',num2str(inh_hz_index));
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        % ignore inhibitions of 20000 and 50000
        
        temp_fields = fields(output);
        c20000 = strfind(fields(output),'c20000');
        c50000 = strfind(fields(output),'c50000');
        for temp_index = 1:length(temp_fields)
            if(isempty(c20000{temp_index})==0 || isempty(c50000{temp_index})==0)
                struct_val = temp_fields{temp_index};
                output = rmfield(output,struct_val);
                
            end
        end
        % transpose data so that strength is by synapse number instead of
        % frequency
        output = transposeSynapticNumberToStrength(output);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        f = figure();
        slope=0;
        for split_index = 1:length(split_fields)
            if (strcmp(split_fields{split_index},strcat('e',num2str(exc_hz_index),'i',num2str(inh_hz_index))))
                top_split_index = split_index;
                split_result = output.(split_fields{split_index});
                [degree_fields, degree_fields_num] = reorderByNum(fields(split_result));
                colors = getColorPalette();
                c=0;
                for degree_index = 1:length(fields(split_result))
                    degree_result = split_result.(degree_fields{degree_index});
                    excitation = degree_result.Excitation;
                    spikes = degree_result.Spikes;
                    top_excitation{degree_index,k} = excitation;
                    top_spikes{degree_index,k} = spikes;

                    if(slope~=0)
                        initial_params = [min(spikes), max(spikes),x50,slope];
                    else
                        initial_params = [];
                    end
                    fixed_params=[min(spikes), spikes(end) , NaN , NaN];
                    % fit to sigmoid
                    [params,stat,x,y] = sigm_fit(log10(excitation),spikes,fixed_params,initial_params,0);
                    x50 = params(3);
                    slope = params(4);
                    % plot result
                    if(toPlot)
                        if(c>length(colors))
                            c=1;
                        else
                            c=c+1;
                        end
                        % fitted sigmoid
                        plots(degree_index)=semilogx(10.^x,y,'Color',colors(c,:));
                        %c = get(plots(degree_index),'Color');
                        hold on;
                        % actual points
                        temp = semilogx(excitation,spikes,'o');
                        set(temp,'Color',colors(c,:));
                    end
                end
                if(toPlot)
                    legend(plots,degree_fields_num);
                end
            end
        end
        split_fields = multiReplaceAll(split_fields,{'e','i'},{'E:',' I:'}); 
    end
    
    % Calculate Normalization Index
    EX = top_excitation{top_split_index,1};
    IN = degree_fields;
    % which is fixed and which is free based on order of keepfunctions
    INxFRFixed = top_spikes(:,1);
    INxFRFree = top_spikes(:,2);
    [Nindex, fig] = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree,...
        'units','# synapses',...
        'format','%.0f',...
        'synapses',1);
    fig = formatFigure(fig);
    f = figure();
    plot(1:length(IN),Nindex);
    ax = gca;
    ax.XTickLabels = multiReplaceAll(IN,{'c'},{''});
    result = Nindex;
end
function result = Data_2016_09_21_17h39_compare_synapses_increasing_synapse_numbers()
    % 10 Hz inhibitory input, increasing excitatory input.
    % increase numbers of synapse to determine gain shift
    % proximal + distal
    % with + without KCC2
    % synapses aren't evenly distributed but clumped into groups (10
    % groups)
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_21_17h39_compare_synapses_increasing_synapse_numbers';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    inh_hz_index = 10;
    exc_hz_index = 10;
    inh_hz_index_str = strcat('d',num2str(inh_hz_index));
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        % transpose data so that strength is by synapse number instead of
        % frequency
        output = transposeSynapticNumberToStrength(output);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        f = figure();
        
        for split_index = 1:length(split_fields)
            if (strcmp(split_fields{split_index},strcat('e',num2str(exc_hz_index),'i',num2str(inh_hz_index))))
                top_split_index = split_index;
                split_result = output.(split_fields{split_index});
                [degree_fields, degree_fields_num] = reorderByNum(fields(split_result));
                for degree_index = 1:length(fields(split_result))
                    degree_result = split_result.(degree_fields{degree_index});
                    excitation = degree_result.Excitation;
                    spikes = degree_result.Spikes;
                    top_excitation{degree_index,k} = excitation;
                    top_spikes{degree_index,k} = spikes;
                    semilogx(excitation,spikes);
                    hold on;
                end
                legend(degree_fields_num);
            end
        end
        split_fields = multiReplaceAll(split_fields,{'e','i'},{'E:',' I:'}); 
    end
    
    % Calculate Normalization Index
    EX = top_excitation{top_split_index,1};
    IN = degree_fields;
    % which is fixed and which is free based on order of keepfunctions
    INxFRFixed = top_spikes(:,1);
    INxFRFree = top_spikes(:,2);
    Nindex = NormalizationIndex(EX,IN,INxFRFixed,INxFRFree);
    result = Nindex;
end
function result = Data_2016_09_15_16h21_10_Hz_increase_inhibitory_synapse_numbers_persistent_excitation()
    % 10 Hz inhibitory input, increasing PERSISTENT excitatory input.
    % increase numbers of synapse to determine gain shift
    % proximal + distal
    % with + without KCC2
    % synapses aren't evenly distributed but clumped into groups (10
    % groups)
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_15_16h21_10_Hz_increase_inhibitory_synapse_numbers_persistent_excitation';
    folderName = foldername;
    % exhaustive set of functions to analyse data
    keepFunctions = {@getKeepDistalNoKCC2};
    SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    result = 1;
end
function result = Data_2016_09_15_15h32_compare_synapses_increasing_synapse_numbers()
    % 10 Hz inhibitory and excitatory input.
    % increase numbers of synapse to determine gain shift
    % proximal + distal
    % with + without KCC2
    % synapses aren't evenly distributed but clumped into groups (10
    % groups)
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_15_15h32_compare_synapses_increasing_synapse_numbers';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    %SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    %CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        %topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        new_output = struct;
        for s = 1:6
            new_split = split_fields{s}(5:end);
            for ns = 1:6
                tmp_struct = output.(split_fields{s+(ns-1)*6}).d10;
                new_output.(new_split).excitation(ns) = str2double(split_fields{s+(ns-1)*6}(2:4));
                new_output.(new_split).spikes(ns) = tmp_struct.Spikes;
            end
        end
        
        output=new_output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f = figure();
        hold on;
        for split_index = 1:length(split_fields)
            split_result = output.(split_fields{split_index});
            excitation = split_result.excitation;
            spikes = split_result.spikes;
            plot(excitation,spikes);
        end
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        result = 1;    
    end
    result = 1;
end
function result = Data_2016_09_09_15h10()
    % 10 Hz inhibitory input, increasing excitatory input.
    % increase numbers of synapse to determine gain shift
    % proximal + distal
    % with + without KCC2
    % synapses aren't evenly distributed but clumped into groups (10
    % groups)
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_09_15h10_10_Hz_increase_inhibitory_synapse_numbers_clumped';
    folderName = foldername;
    % exhaustive set of functions to analyse data
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2,@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    %CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    %keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    %CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    %keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    %CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    %keepFunctions = {@getKeepDistal,@getKeepproximal};
    %result = splitOnKCC2DegreeOfName(foldername,keepFunctions);
    result = 1;
end
function result = Data_2016_09_09_compare_base_steady_cli()
    % run distal pattern of inhibition for matched input, but using
    % distal's optimal synaptic numbers
    % compare using cli of 5 (set) versus ~4.25 (steady-state)
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_09_13h51_compare_synapses_distal_balanced_synapses';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalNoKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    keepnames{1} = 'steadystate';
    funcnames{1} = 'steadystate';
    keepnames{2} = 'set';
    funcnames{2} = 'set';
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        if(k==2)
            % compare to files in another folder
            foldername = 'Hz/Data_2016_09_09_10h25_compare_synapsesdistal_balanced_synapses';
        else
            foldername = folderName;
        end
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(funcnames{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(funcnames{k}).excitation = excitation;
        topoutput.(funcnames{k}).spikes = spikes;
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(funcnames{k}).excitation(split_index,:),topoutput.(funcnames{k}).spikes(split_index,:))
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(funcnames{1}).excitation(split_index,:),topoutput.(funcnames{1}).spikes(split_index,:),'-','Color',colors(split_index,:));
        plot(topoutput.(funcnames{2}).excitation(split_index,:),topoutput.(funcnames{2}).spikes(split_index,:),'--','Color',colors(split_index,:));
    end
    plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    title('Distal (-: steady-state cli, --: set cli)');
    saveFigs('folderName',{'with and without set cli'},[f3]);
    
end
function result = Data_2016_09_09_13h51_compare_synapses_distal_balanced_synapses()
    % run distal pattern of inhibition for matched input, but using
    % distal's optimal synaptic numbers
    % difference from below method is equal steady-state chloride for all
    % runs
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_09_13h51_compare_synapses_distal_balanced_synapses';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        if(strcmp(funcnames{k},'getKeepDistalKCC2'))
            % compare to files in another folder
            foldername = 'Hz/Data_2016_09_09_10h25_compare_synapsesdistal_balanced_synapses';
        else
            foldername = folderName;
        end
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(funcnames{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(funcnames{k}).excitation = excitation;
        topoutput.(funcnames{k}).spikes = spikes;
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(funcnames{k}).excitation(split_index,:),topoutput.(funcnames{k}).spikes(split_index,:))
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(funcnames{1}).excitation(split_index,:),topoutput.(funcnames{1}).spikes(split_index,:),'-','Color',colors(split_index,:));
        plot(topoutput.(funcnames{2}).excitation(split_index,:),topoutput.(funcnames{2}).spikes(split_index,:),'--','Color',colors(split_index,:));
    end
    plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    title('Distal (-: no chloride dynamics, --: with chloride dynamics)');
    saveFigs('folderName',{'with and without KCC2'},[f3]);
    
end
function result = Data_2016_09_09_10h25_compare_synapsesdistal_balanced_synapses()
    % run distal pattern of inhibition for matched input, but using
    % distal's optimal synaptic numbers
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_09_10h25_compare_synapsesdistal_balanced_synapses';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}, funcnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(funcnames{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(funcnames{k}).excitation = excitation;
        topoutput.(funcnames{k}).spikes = spikes;
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(funcnames{k}).excitation(split_index,:),topoutput.(funcnames{k}).spikes(split_index,:))
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(funcnames{1}).excitation(split_index,:),topoutput.(funcnames{1}).spikes(split_index,:),'-','Color',colors(split_index,:));
        plot(topoutput.(funcnames{2}).excitation(split_index,:),topoutput.(funcnames{2}).spikes(split_index,:),'--','Color',colors(split_index,:));
    end
    plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    title('Distal (-: no chloride dynamics, --: with chloride dynamics)');
    saveFigs('folderName',{'with and without KCC2'},[f3]);
    
end
function result = Data_2016_09_07_increase_inhibitory_synapses_numbers()
    % 10 Hz inhibitory input, increasing excitatory input.
    % increase numbers of synapse to determine gain shift
    % proximal + distal
    % with + without KCC2
    global folderName
    whos
    foldername = 'Hz/Data_2016_09_07_increase_inhibitory_synapses_numbers';
    folderName = foldername;
    % exhaustive set of functions to analyse data
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2,@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    SplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    CompareKeepsSplitOnNameDegreeOfInhibition(foldername,keepFunctions);
    keepFunctions = {@getKeepDistal,@getKeepproximal};
    result = splitOnKCC2DegreeOfName(foldername,keepFunctions);
    
end
function result = Data_2016_08_check_pop_out_control_run()
    % Run with matched input at different numbers of synapses with and
    % without KCC2 for distal and proximal with control curves (inhibition
    % at 0 Hz) to see expected excitation IO curve.
    global folderName
    whos
    %% Proximal versions
    foldername = 'Hz/Data_2016_08_13_14h19_control_run';
    folderName = foldername;
    savename = 'Hz/Data_2016_08_check_pop_out_control_run';
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    keepStrings = {'getKeepproximalNoKCC2','getKeepproximalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        %topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitationControl{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikesControl{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           %plot(excitation{split_index,:},spikes{split_index,:})
        end
        topoutput.(keepStrings{k}).excitationControl = excitationControl;
        topoutput.(keepStrings{k}).spikesControl = spikesControl;
        %plot(0:5:5*length(excitationControl{1,:})-1,0:5:5*length(spikesControl{1,:})-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        result = 1;    
    end
    
    
    % Run with matched input at different numbers of synapses with and
    % without KCC2
    foldername = 'Hz/Data_2016_08_06_19h01_check_pop_out';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    keepStrings = {'getKeepproximalNoKCC2','getKeepproximalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        %topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
%         f(k) = figure();
%         hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           %plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(keepStrings{k}).excitation = excitation;
        topoutput.(keepStrings{k}).spikes = spikes;
        %plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        result = 1;    
    end
    %saveFigs('folderName',keepnames,f);
    keepnamesControl = strcat(keepnames,' control');
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(keepStrings{k}).excitation(split_index,:),topoutput.(keepStrings{k}).spikes(split_index,:))
           plot(topoutput.(keepStrings{k}).excitationControl{split_index,:},topoutput.(keepStrings{k}).spikesControl{split_index,:})
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend([keepnames(1), keepnamesControl(1), keepnames(2), keepnamesControl(2)]);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs(savename,strcat('proximal_',split_fields),f2);
    
    
    %% Distal versions
    foldername = 'Hz/Data_2016_08_13_14h19_control_run';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    keepStrings = {'getKeepDistalNoKCC2','getKeepDistalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        %topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitationControl{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikesControl{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           %plot(excitation{split_index,:},spikes{split_index,:})
        end
        topoutput.(keepStrings{k}).excitationControl = excitationControl;
        topoutput.(keepStrings{k}).spikesControl = spikesControl;
        %plot(0:5:5*length(excitationControl{1,:})-1,0:5:5*length(spikesControl{1,:})-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        result = 1;    
    end
    
    
    % Run with matched input at different numbers of synapses with and
    % without KCC2
    foldername = 'Hz/Data_2016_08_05_16h30_check_pop_out';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    keepStrings = {'getKeepDistalNoKCC2','getKeepDistalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        %topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
%         f(k) = figure();
%         hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           %plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(keepStrings{k}).excitation = excitation;
        topoutput.(keepStrings{k}).spikes = spikes;
        %plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        result = 1;    
    end
    %saveFigs('folderName',keepnames,f);
    keepnamesControl = strcat(keepnames,' control');
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(keepStrings{k}).excitation(split_index,:),topoutput.(keepStrings{k}).spikes(split_index,:))
           plot(topoutput.(keepStrings{k}).excitationControl{split_index,:},topoutput.(keepStrings{k}).spikesControl{split_index,:})
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend([keepnames(1), keepnamesControl(1), keepnames(2), keepnamesControl(2)]);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs(savename,strcat('distal_',split_fields),f2);
    
    
end
function result = Data_2016_08_13_14h19_control_run()
    % Run with matched input at different numbers of synapses with and
    % without KCC2
    global folderName
    whos
    foldername = 'Hz/Data_2016_08_13_14h19_control_run';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2,@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    keepStrings = {'getKeepproximalNoKCC2','getKeepproximalKCC2','getKeepDistalNoKCC2','getKeepDistalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation{split_index,:},spikes{split_index,:})
        end
        topoutput.(keepStrings{k}).excitation = excitation;
        topoutput.(keepStrings{k}).spikes = spikes;
        plot(0:5:5*length(excitation{1,:})-1,0:5:5*length(spikes{1,:})-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(keepStrings{k}).excitation{split_index,:},topoutput.(keepStrings{k}).spikes{split_index,:})
        end
        plot(0:5:5*length(excitation{1,:})-1,0:5:5*length(spikes{1,:})-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
end
function result = Data_2016_08_06_19h01_check_pop_out()
    % Run with matched input at different numbers of synapses with and
    % without KCC2
    global folderName
    whos
    foldername = 'Hz/Data_2016_08_06_19h01_check_pop_out';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2,@getKeepproximalKCC2};
    keepStrings = {'getKeepproximalNoKCC2','getKeepproximalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(keepStrings{k}).excitation = excitation;
        topoutput.(keepStrings{k}).spikes = spikes;
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    % plot separate figures comparing with and without KCC2
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(keepStrings{k}).excitation(split_index,:),topoutput.(keepStrings{k}).spikes(split_index,:))
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing woth and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(keepStrings{1}).excitation(split_index,:),topoutput.(keepStrings{1}).spikes(split_index,:),'-','Color',colors(split_index,:));
        plot(topoutput.(keepStrings{2}).excitation(split_index,:),topoutput.(keepStrings{2}).spikes(split_index,:),'--','Color',colors(split_index,:));
    end
    plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    title('Proximal (-: no chloride dynamics, --: with chloride dynamics)');
    saveFigs('folderName',{'proximal with and without KCC2'},[f3]);
end
function result = Data_2016_08_05_16h30_check_pop_out()
    % Run with matched input at different numbers of synapses with and
    % without KCC2
    global folderName
    whos
    foldername = 'Hz/Data_2016_08_05_16h30_check_pop_out';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2,@getKeepDistalKCC2};
    keepStrings = {'getKeepDistalNoKCC2','getKeepDistalKCC2'};
    [splitReg,splitName,trimBeginSplit,trimEndSplit]=getSplitOnNameInPt();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        topoutput.(keepStrings{k}) = output;
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        f(k) = figure();
        hold on;
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA(split_index,degree_index) = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
           plot(excitation(split_index,:),spikes(split_index,:))
        end
        topoutput.(keepStrings{k}).excitation = excitation;
        topoutput.(keepStrings{k}).spikes = spikes;
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        legend(split_fields);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(keepname);
        result = 1;    
    end
    saveFigs('folderName',keepnames,f);
    for split_index = 1:length(split_fields)
        f2(split_index) = figure();
        hold on;
        for k=1:length(keeps)
           plot(topoutput.(keepStrings{k}).excitation(split_index,:),topoutput.(keepStrings{k}).spikes(split_index,:))
        end
        plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
        legend(keepnames);
        xlabel('Input (Hz)');
        ylabel('Output (Hz)');
        title(split_fields{split_index});
    end
    saveFigs('folderName',split_fields,f2);
    % plot one figure comparing with and without KCC2
    f3 = figure();
    hold on;
    colors = get(groot,'DefaultAxesColorOrder');
    for split_index = 1:length(split_fields)
        p(split_index)=plot(topoutput.(keepStrings{1}).excitation(split_index,:),topoutput.(keepStrings{1}).spikes(split_index,:),'-','Color',colors(split_index,:));
        plot(topoutput.(keepStrings{2}).excitation(split_index,:),topoutput.(keepStrings{2}).spikes(split_index,:),'--','Color',colors(split_index,:));
    end
    plot(0:5:5*length(excitation(1,:))-1,0:5:5*length(spikes(1,:))-1,'ko--');
    legend(p,split_fields);
    xlabel('Input (Hz)');
    ylabel('Output (Hz)');
    title('Distal (-: no chloride dynamics, --: with chloride dynamics)');
    saveFigs('folderName',{'with and without KCC2'},[f3]);
end
function result = Data_2016_08_02_11h27_run_synapse_types()
    % Run with synaptic type of 3 (tonic excitation, synaptic inhibition),
    % but also see what an increase in synaptic weight does
    global folderName
    whos
    foldername = 'Hz/Data_2016_08_02_11h27_run_synapse_types';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    [split,splitName,trimBeginSplit,trimEndSplit]=getSplitOnWeights();
    [degree,degreeName,trimBeginDegree,trimEndDegree]=getDegreeOfInhibition();                  
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',trimBeginSplit,'trimEndSplit',trimEndSplit,'trimBeginDegree',trimBeginDegree,'trimEndDegree',trimEndDegree);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields_sub)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
           end
           degree_fields{split_index} = degree_fields_sub;
           degree_fields_num{split_index} = degree_fields_sub_num;
        end
        split_fields = multiReplaceAll(split_fields,{'s','c'},{'E:',' I:'});
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split','rowtitle',...
                splitName,'columntitle',degreeName);
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree','rowtitle',...
                splitName,'columntitle',degreeName);
        result = 1;    
    end
end
function result = Data_2016_07_29_18h49_run_synapse_types()
    % Run with 10 synapses of each type, but with different types of
    % synapses.
    %   0: both synaptic
    %   1: both tonic
    %   2: synaptic excitation, tonic inhibition
    %   3: tonic excitation, synaptic inhibition
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_29_18h49_run_synapse_types';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{k};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',1,'trimEndSplit',5,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        allowedCell={};
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields_sub,degree_fields_sub_num] = reorderByNum(fields(split_result),cutEnd);
           %# When a mismatch between degree fields occurs, need to reduce to common values
           if(isempty(allowedCell))
              allowedCell = degree_fields_sub;
           end
           used_index=1;
           for degree_index = 1:length(degree_fields_sub)
               allowed=false;
               for ai = 1:length(allowedCell)
                  if(strcmp(degree_fields_sub{degree_index},allowedCell{ai}))
                     allowed=true;
                     break;
                  end
               end
               if(allowed)
                   excitation{split_index,used_index} = split_result.(degree_fields_sub{degree_index}).Excitation; 
                   spikes{split_index,used_index} = split_result.(degree_fields_sub{degree_index}).Spikes;
                   somaEGABA{split_index,used_index} = split_result.(degree_fields_sub{degree_index}).Soma_EGABA;
                   ldendEGABA{split_index,used_index} = split_result.(degree_fields_sub{degree_index}).ldend_EGABA;
                   bdendEGABA{split_index,used_index} = split_result.(degree_fields_sub{degree_index}).bdend_EGABA_;
                   usedDF{used_index}=degree_fields_sub{used_index};
                   usedDFN{used_index}=degree_fields_sub_num{used_index};
                   used_index=used_index+1;
               end
           end
           degree_fields{split_index} = usedDF;
           degree_fields_num{split_index} = usedDFN;
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_29_12h17_Extended_exc_range()
    % Run previous (120:10 synapses) over broader range of excitation
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_29_12h17_Extended_exc_range';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',1,'trimEndSplit',5,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_28_14h01_run_synapse_types()
    % Run with 10 synapses of each type, but with tonic inhibition and 
    % synaptic excitation
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_28_14h01_run_synapse_types';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',1,'trimEndSplit',5,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_28_10h18_run_num_synapses()
    % Run with 10 synapses of each type, but with "best" weights (see
    % default model)
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_28_10h18_run_num_synapses';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',2,'trimEndSplit',6,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_27_12h12()
    % Same number of synapses, different weights
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_27_12h12';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',2,'trimEndSplit',6,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_26_10h28pCl()
    % Compare pCl values
    % Non KCC2 only, for proof of principle
    global folderName    
    whos
    foldername = 'Hz/Data_2016_07_26_10h28(pCl)';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',1,'trimEndSplit',5,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        CompareInnerNormalizationIndices(excitation,spikes,split_fields_num,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields_num,...
                degree_fields_num,keepname,'orderby','degree');
        result = 1;    
    end
end
function result = Data_2016_07_25_11h19()
    % Compare synapse number behaviours at higher levels of excitation
    % (DISTAL)
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_25_11h19';
    folderName = foldername;
    keepFunctions = {@getKeepDistalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on synapse numbers
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',2,'trimEndSplit',6,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
        
        split_fields = multiReplace(split_fields,{'s','c'},{'E:',' I:'});
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','column');

        result = 1;    
    end
end
function result = Data_2016_07_22_10h48()
    % Compare synapse number behaviours at higher levels of excitation
    % (PROXIMAL)
    % Non KCC2 only, for proof of principle
    global folderName
    whos
    foldername = 'Hz/Data_2016_07_22_10h48';
    folderName = foldername;
    keepFunctions = {@getKeepproximalNoKCC2};
    for kf = 1:length(keepFunctions)
        [keeps{kf}, keepnames{kf}] = keepFunctions{kf}();
    end
    for k = 1:length(keeps)
        keep = keeps{k};
        keepname = keepnames{1};
        split='_.+_InPt';                   %split on values between name and inhibition
        degree='InPt\d+_';                  %degree of inhibition strength
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',2,'trimEndSplit',6,'trimBeginDegree',4,'trimEndDegree',1);
        [split_fields, split_fields_num] = reorderByNum(fields(output));
        
        for split_index = 1:length(split_fields)
           split_result = output.(split_fields{split_index});
           cutEnd=0;    
           [degree_fields,degree_fields_num] = reorderByNum(fields(split_result),cutEnd);
           for degree_index = 1:length(degree_fields)
               %columns: same degree, different split (:,df)
               %rows: same split, different degree (rf,:)
               excitation{split_index,degree_index} = split_result.(degree_fields{degree_index}).Excitation; 
               spikes{split_index,degree_index} = split_result.(degree_fields{degree_index}).Spikes;
               somaEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).Soma_EGABA;
               ldendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).ldend_EGABA;
               bdendEGABA{split_index,degree_index} = split_result.(degree_fields{degree_index}).bdend_EGABA_;
           end
        end
       CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','split');
        CompareInnerNormalizationIndices(excitation,spikes,split_fields,...
                degree_fields_num,keepname,'orderby','degree');
        
        
        result = 1;    
    end
end
function result = Data_2016_07_19_15h31()
    %Compare slight diameter changes (to extend still)
    whos
    
   %% diam comparison
    foldername = 'Hz/Data_2016_07_19_15h31';
    keeps = {getKeepproximalKCC2(),getKeepproximalNoKCC2()};
    for k = 1:length(keeps)
        keep = keeps{k};
        split='InPt\d+_';                       %split on inhibition strength
        degree='_[.\d]+.*InPt';                 %degrees of diam
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',4,'trimEndSplit',1,'trimEndDegree',4);
        [inhib_fields, inhib_fields_num] = reorderByNum(fields(output));
        
        
        for inhf = 1:length(inhib_fields)
           inhib_result = output.(inhib_fields{inhf});
           cutEnd=0;    
           [diam_fields,diam_fields_num] = reorderByNum(fields(inhib_result),cutEnd);
           for diam = 1:length(diam_fields)
               %rows: same inhibition, different diam (rf,:)
               %columns: same KCC2Pa, different inhibition (:,df)
               excitation{inhf,diam} = inhib_result.(diam_fields{diam}).Excitation;
               spikes{inhf,diam} = inhib_result.(diam_fields{diam}).Spikes;
               somaEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).Soma_EGABA;
               ldendEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).ldend_EGABA;
               bdendEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).bdend_EGABA_;
           end
        end
        %%TODO: put fixed cl in matrices for nindex calc
         for inhf = 1:length(inhib_fields)
             Nindex2{inhf} = InnerNormalizationIndex(excitation{inhf},spikes(inhf,:),strcat(regexprep(keep,'_','-'),inhib_fields{inhf}));
    %         %Nindex_soma_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},somaEGABA(inhf,:),inhib_fields{inhf});
    %         %Nindex_ldend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},ldendEGABA(inhf,:),inhib_fields{inhf});
    %         %Nindex_bdend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},bdendEGABA(inhf,:),inhib_fields{inhf});
         end
        for diam = 1:length(diam_fields)
            %Nindex between inhibition strengths for each diam
            Nindex{diam} = InnerNormalizationIndex(excitation{diam},spikes(:,diam),strcat(regexprep(keep,'_','-'),diam_fields{diam}));
        end
        figure;
        hold on;
        for n = 1:length(Nindex)
            plotN=Nindex{n}(~isnan(Nindex{n}));
            %plotNSE=Nindex_soma_EGBABA{n}(~isnan(Nindex_soma_EGBABA{n}));
            %plotNLE=Nindex_ldend_EGBABA{n}(~isnan(Nindex_ldend_EGBABA{n}));
            %plotNBE=Nindex_bdend_EGBABA{n}(~isnan(Nindex_bdend_EGBABA{n}));
            %plotNE = bdendEGABA(inhf,:);
            %[plotN,plotNE] = matchNindices(plotN,plotNE);
            %if(length(plotN)==length(inhib_fields))
                plot(1:length(plotN),plotN,'-+');
                %legplot = plot(plotNE,plotN,'-+');
                %plot(1:length(plotNSE),plotNSE,'-or');
                %plot(1:length(plotNLE),plotNLE,'-xg');
                %plot(1:length(plotNSE),plotNSE,'*');
            %end
        end
        legend(diam_fields_num);
        set(gca,'XTick',linspace(1,length(inhib_fields),length(inhib_fields)));
        set(gca,'XTickLabel',inhib_fields_num);
        xlabel('inhibition strength (x baseline conductance)');
        ylabel('"Normalisation" index');
        if(keep == getKeepproximalNoKCC2())
            title(multiReplace(keep,{'[^K][^C]+[^2][^_].*'},{'NO KCC2'}));
        end
        
        
        result = Nindex;    
    end
end
function result = Data_2016_07_18_14h56()
    %% use of GABA and NMDA.mod
    %  higher values of frequencies used, check response direction
    % noise     1       1
    % # syns:   50     140
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-18_14h56';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, distal, distal_KCC2, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    % DISTAL
    distal_fields = fields(distal);
    EX_distal = distal.(distal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(distal_fields)
        INx3{x} = distal.(distal_fields{x}).Spikes;
        INx4{x} = distal_KCC2.(distal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_distal] = NormalizationIndex(EX_distal,distal_fields,INx3,INx4);
    disp(Nindex_distal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,distal_fields,EX_distal,INx3,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,distal_fields,EX_distal,INx4,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_distal;
end
function result = Data_2016_07_18_11h58()
    %% use of GABA and NMDA.mod
    %  higher values of frequencies used, check response direction
    % noise     1       1
    % # syns:   10     120
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-18_11h58';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, distal, distal_KCC2, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    % DISTAL
    distal_fields = fields(distal);
    EX_distal = distal.(distal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(distal_fields)
        INx3{x} = distal.(distal_fields{x}).Spikes;
        INx4{x} = distal_KCC2.(distal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_distal] = NormalizationIndex(EX_distal,distal_fields,INx3,INx4);
    disp(Nindex_distal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,distal_fields,EX_distal,INx3,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,distal_fields,EX_distal,INx4,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_distal;
end
function result = Data_2016_07_17_22h59()
    %% use of GABA and NMDA.mod
    %  higher values of frequencies used, check response direction
    % noise     0       0.1
    % # syns:  130    180
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-17_22h59';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, ~, ~, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_proximal;
end
function result = Data_2016_07_17_11h10()
    %% use of GABA and NMDA.mod
    %  higher values of frequencies used, check response direction
    % noise     0       0
    % # syns:  130    180
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-17_11h10';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, ~, ~, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_proximal;
end
function result = Data_2016_07_14_10h45()
    %% use of GABA and NMDA.mod
    %  higher values of frequencies used, check response direction
    % noise     1       1
    % # syns:  130    180
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-14_10h45';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, ~, ~, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_proximal;
end
function result = Data_2016_07_05_13h06()
    %% finer lower inhibition and excitation values
    % noise with 50 pS weight
    global toPlot toSave
    whos
    %% KCC2 Comparison
    foldername = 'Hz/Data-2016-07-05_13h06';
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;

    [proximal, proximal_KCC2, ~, ~, ~, ~]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end

    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    disp(Nindex_proximal);
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    desiredFR = 5;
    units='Hz';
    titleName='Proximal no Cl dynamics';
    [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
    titleName='Proximal with Cl dynamics';
    [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
    [fig1,fig2] = adjustLims(fig1,fig2);
    h(1) = fig1; h(2) = fig2;
    newfig = combineFigures(h);
    if(save)
        saveNames = {'followCurveNoCl','followCurveWithCl'};
        saveFigs(foldername,saveNames,[fig1,fig2]);
    end
    result = Nindex_proximal;
end
function result = Data_2016_06_04_14h20()
    %% finer lower inhibition and excitation values
    % noise with 50 weight
    global toPlot toSave
    whos
    old = cd('Hz');
    %% KCC2 Comparison
    foldername = 'Data-2016-06-04_14h20';
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    extras = {'_IW50'};
    for extraCell = extras
        extra = extraCell{1};
        [proximal, proximal_KCC2, ~, ~, ~, ~]...
            = compareSynapticPatternAndKCC2(foldername,extra,...
                                            prox,dist,inter,toplot,save);

        % Get excitation and spikes (at axon)
        proximal_fields = fields(proximal);
        EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
        for x=1:length(proximal_fields)
            INx1{x} = proximal.(proximal_fields{x}).Spikes;
            INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
        end

        % Generate Nindex
        [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
        disp(Nindex_proximal);
        % Check 'desired firing rate' curve
        inhibitonX = 1;
        desiredFR = 5;
        units='Hz';
        titleName=strcat('Proximal no Cl dynamics',extra);
        [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
        titleName=strcat('Proximal with Cl dynamics',extra);
        [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
        [fig1,fig2] = adjustLims(fig1,fig2);
        h(1) = fig1; h(2) = fig2;
        newfig = combineFigures(h);
        if(save)
            saveNames = strcat({'followCurveNoCl','followCurveWithCl'},extra);
            saveFigs(foldername,saveNames,[fig1,fig2]);
        end
        result = Nindex_proximal;
    end
    cd(old);
end
function result = Data_2016_06_03_19h14()
    %% finer lower inhibition and excitation values
    % No noise for netstim firing
    global toPlot toSave
    whos
    old = cd('Hz');
    %% KCC2 Comparison
    foldername = 'Data-2016-06-03_19h14';
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    cd(old);
    extras = {'_IW10','_IW2'};
    for extraCell = extras
        extra = extraCell{1};
        [proximal, proximal_KCC2, ~, ~, ~, ~]...
            = compareSynapticPatternAndKCC2(foldername,extra,...
                                            prox,dist,inter,toplot,save);

        % Get excitation and spikes (at axon)
        proximal_fields = fields(proximal);
        EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
        for x=1:length(proximal_fields)
            INx1{x} = proximal.(proximal_fields{x}).Spikes;
            INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
        end

        % Generate Nindex
        [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
        disp(Nindex_proximal);
        % Check 'desired firing rate' curve
        inhibitonX = 1;
        desiredFR = 5;
        units='Hz';
        titleName=strcat('Proximal no Cl dynamics',extra);
        [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
        titleName=strcat('Proximal with Cl dynamics',extra);
        [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
        h(1) = fig1; h(2) = fig2;
        newfig = combineFigures(h);
        if(save)
            saveNames = strcat({'followCurveNoCl','followCurveWithCl'},extra);
            saveFigs(foldername,saveNames,[fig1,fig2]);
        end
        result = Nindex_proximal;
        
    end
    cd(old);
end
function result = Data_2016_06_03_14h48()
    %% finer lower inhibition and excitation values
    global toPlot toSave
    whos
    cd('Hz');
    %% KCC2 Comparison
    foldername = 'Data-2016-06-03_14h48';
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    extra='_IW10'; %# _IW10 or _IW2
    extras = {'_IW10','_IW2'};
    for extraCell = extras
        extra = extraCell{1};
        [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
            = compareSynapticPatternAndKCC2(foldername,extra,...
                                            prox,dist,inter,toplot,save);

        % Get excitation and spikes (at axon)
        proximal_fields = fields(proximal);
        EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
        for x=1:length(proximal_fields)
            INx1{x} = proximal.(proximal_fields{x}).Spikes;
            INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
        end

        % Generate Nindex
        [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
        disp(Nindex_proximal);
        % Check 'desired firing rate' curve
        inhibitonX = 1;
        desiredFR = 6;
        units='Hz';
        titleName=strcat('Proximal no Cl dynamics',extra);
        [result, fig1] = followDesired(desiredFR,proximal_fields,EX_proximal,INx1,inhibitonX,titleName,units);
        titleName=strcat('Proximal with Cl dynamics',extra);
        [result, fig2] = followDesired(desiredFR,proximal_fields,EX_proximal,INx2,inhibitonX,titleName,units);
        if(save)
            saveNames = strcat({'followCurveNoCl','followCurveWithCl'},extra);
            saveFigs(foldername,saveNames,[fig1,fig2]);
        end
        result = Nindex_proximal;
    end
end
function result = Data_2016_06_02_16h12()
    %% finer lower inhibition and excitation values
    % The inhibition weights are too weak
    global toPlot toSave
    whos
    cd('Hz');
    %% KCC2 Comparison
    foldername = 'Data-2016-06-02_16h12';
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        distal_fields = fields(distal);
        EX_distal = distal.(distal_fields{1}).Excitation;           %Excitation values are the same for all distal runs
    end
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end
    
    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    [result, fig1] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
    [result, fig2] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    if(save)
       saveFigs(foldername,{'graph_followCurveNoCl','graph_followCurveWithCl'},[fig1,fig2])
    end
    result = Nindex_proximal;
end
function result = Data_2016_06_02_15h00()
    %% finer lower inhibition and excitation values
    % The inhibition weights are too weak
    global toPlot toSave
    whos
    cd('Hz');
    %% KCC2 Comparison
    foldername = 'Data-2016-06-02_15h00';
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2(foldername,'',...
                                        prox,dist,inter,toplot,save);
    
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        distal_fields = fields(distal);
        EX_distal = distal.(distal_fields{1}).Excitation;           %Excitation values are the same for all distal runs
    end
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end
    
    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    [result, fig1] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
    [result, fig2] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    if(save)
       saveFigs(foldername,{'graph_followCurveNoCl','graph_followCurveWithCl'},[fig1,fig2])
    end
    result = Nindex_proximal;
end
function result = Data_2016_05_29_moreInhib()
    global toPlot toSave
    whos
    old = cd('Hz');
    %% KCC2 Comparison
    prox = 1;
    dist=0;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2('Data-2016-05-29(moreInhib)','1)',...
                                        prox,dist,inter,toplot,save);
    cd(old);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        distal_fields = fields(distal);
        EX_distal = distal.(distal_fields{1}).Excitation;           %Excitation values are the same for all distal runs
    end
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end
    
    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    
    if(dist) 
        x_index=1;
        for x=1:length(distal_fields)
            INx3{x} = distal.(distal_fields{x}).Spikes;
            INx4{x} = distal_KCC2.(distal_fields{x}).Spikes;
            if(not(and(distal_KCC2.(distal_fields{x}).Excitation(1)==0,distal_KCC2.(distal_fields{x}).Spikes(1)>5)))
                % Filter out excitatory inhibition
                INx3_filtered{x_index} = distal.(distal_fields{x}).Spikes;
                INx4_filtered{x_index} = distal_KCC2.(distal_fields{x}).Spikes;
                x_index=x_index+1;
            end
        end
       [Nindex_distal] = NormalizationIndex(EX_distal,distal_fields,INx3,INx4);  
    end 
    result = Nindex_proximal;
end
function result = Data_2016_05_25_pCl0_8Hz()
    global toPlot toSave
    whos
    old = cd('Hz');
    %% KCC2 Comparison
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2('Data-2016-05-25(pCl0.8)','5)',...
                                        prox,dist,inter,toplot,save);
    cd(old);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        distal_fields = fields(distal);
        EX_distal = distal.(distal_fields{1}).Excitation;           %Excitation values are the same for all distal runs
    end
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end
    
    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    
    if(dist) 
        x_index=1;
        for x=1:length(distal_fields)
            INx3{x} = distal.(distal_fields{x}).Spikes;
            INx4{x} = distal_KCC2.(distal_fields{x}).Spikes;
            if(not(and(distal_KCC2.(distal_fields{x}).Excitation(1)==0,distal_KCC2.(distal_fields{x}).Spikes(1)>5)))
                % Filter out excitatory inhibition
                INx3_filtered{x_index} = distal.(distal_fields{x}).Spikes;
                INx4_filtered{x_index} = distal_KCC2.(distal_fields{x}).Spikes;
                x_index=x_index+1;
            end
        end
       [Nindex_distal] = NormalizationIndex(EX_distal,distal_fields,INx3,INx4);  
       %[result, fig] = followDesired(5,distal_fields,EX_distal,INx3,inhibitonX);
        %[result, fig] = followDesired(5,distal_fields,EX_distal,INx4,inhibitonX);
    end 
    result = Nindex_proximal;
end
function result = Data_2016_05_25_pCl0_8()
    global toPlot toSave
    whos
    old = cd('Andy+Chris');
    %% KCC2 Comparison
    prox = 1;
    dist=1;
    inter = 0;
    toplot = toPlot;
    save=toSave;
    [proximal, proximal_KCC2, distal, distal_KCC2, internrn, internrn_KCC2]...
        = compareSynapticPatternAndKCC2('Data-2016-05-25(pCl0.8)','',...
                                        prox,dist,inter,toplot,save);
    cd(old);
    % Get excitation and spikes (at axon)
    proximal_fields = fields(proximal);
    EX_proximal = proximal.(proximal_fields{1}).Excitation;     %Excitation values are the same for all proximal runs
    if(dist)
        distal_fields = fields(distal);
        EX_distal = distal.(distal_fields{1}).Excitation;           %Excitation values are the same for all distal runs
    end
    for x=1:length(proximal_fields)
        INx1{x} = proximal.(proximal_fields{x}).Spikes;
        INx2{x} = proximal_KCC2.(proximal_fields{x}).Spikes;
    end
    
    % Generate Nindex
    [Nindex_proximal] = NormalizationIndex(EX_proximal,proximal_fields,INx1,INx2);
    
    % Check 'desired firing rate' curve
    inhibitonX = 1;
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx1,inhibitonX);
    [result, fig] = followDesired(5,proximal_fields,EX_proximal,INx2,inhibitonX);
    
    if(dist) 
        x_index=1;
        for x=1:length(distal_fields)
            INx3{x} = distal.(distal_fields{x}).Spikes;
            INx4{x} = distal_KCC2.(distal_fields{x}).Spikes;
            if(not(and(distal_KCC2.(distal_fields{x}).Excitation(1)==0,distal_KCC2.(distal_fields{x}).Spikes(1)>5)))
                % Filter out excitatory inhibition
                INx3_filtered{x_index} = distal.(distal_fields{x}).Spikes;
                INx4_filtered{x_index} = distal_KCC2.(distal_fields{x}).Spikes;
                x_index=x_index+1;
            end
        end
       [Nindex_distal] = NormalizationIndex(EX_distal,distal_fields,INx3,INx4);  
       %[result, fig] = followDesired(5,distal_fields,EX_distal,INx3,inhibitonX);
        %[result, fig] = followDesired(5,distal_fields,EX_distal,INx4,inhibitonX);
    end 
    result = Nindex_proximal;
end
function result = Data_2016_05_18_KCC2Pa()
    old=cd('Andy+Chris');
    whos
    %% KCC2Pa comparison
    foldername = 'Data-2016-05-18(KCC2Pa)';
    keep='proximal*';                      % proximal
    split='InPt[.\d]+(';                    %split on inhibition strength
    degree='_[.\d]+.*InPt';                 %degrees of diam
    output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',4,'trimEndSplit',1,'trimEndDegree',4);
    cd(old);
    [inhib_fields, inhib_fields_num] = reorderByNum(fields(output));
    
    for inhf = 1:length(inhib_fields)
       inhib_result = output.(inhib_fields{inhf});
       cutEnd=2;    % cut the 'e5' at end
       [kcc2pa_fields,kcc2pa_fields_num] = reorderByNum(fields(inhib_result),cutEnd);
       for kccf = 1:length(kcc2pa_fields)
           %rows: same inhibition, different KCC2Pa (rf,:)
           %columns: same KCC2Pa, different inhibition (:,df)
           excitation{inhf,kccf} = inhib_result.(kcc2pa_fields{kccf}).Excitation;
           spikes{inhf,kccf} = inhib_result.(kcc2pa_fields{kccf}).Spikes;
           somaEGABA{inhf,kccf} = inhib_result.(kcc2pa_fields{kccf}).Soma_EGABA;
           ldendEGABA{inhf,kccf} = inhib_result.(kcc2pa_fields{kccf}).ldend_EGABA;
           bdendEGABA{inhf,kccf} = inhib_result.(kcc2pa_fields{kccf}).bdend_EGABA_;
       end
    end
    %%TODO: put fixed cl in matrices for nindex calc
%     for inhf = 1:length(inhib_fields)
%         Nindex{inhf} = InnerNormalizationIndex(excitation{inhf},spikes(inhf,:),inhib_fields{inhf});
%         %Nindex_soma_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},somaEGABA(inhf,:),inhib_fields{inhf});
%         %Nindex_ldend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},ldendEGABA(inhf,:),inhib_fields{inhf});
%         %Nindex_bdend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},bdendEGABA(inhf,:),inhib_fields{inhf});
%     end
    for kccf = 1:length(kcc2pa_fields)
        %Nindex between inhibition strengths for each KCC2Pa strength
        Nindex{kccf} = InnerNormalizationIndex(excitation{kccf},spikes(:,kccf),kcc2pa_fields{kccf});
    end
    figure;
    hold on;
    for n = 1:length(Nindex)
        plotN=Nindex{n}(~isnan(Nindex{n}));
        %plotNSE=Nindex_soma_EGBABA{n}(~isnan(Nindex_soma_EGBABA{n}));
        %plotNLE=Nindex_ldend_EGBABA{n}(~isnan(Nindex_ldend_EGBABA{n}));
        %plotNBE=Nindex_bdend_EGBABA{n}(~isnan(Nindex_bdend_EGBABA{n}));
        %plotNE = bdendEGABA(inhf,:);
        %[plotN,plotNE] = matchNindices(plotN,plotNE);
        %if(length(plotN)==length(inhib_fields))
            plot(1:length(plotN),plotN,'-+');
            %legplot = plot(plotNE,plotN,'-+');
            %plot(1:length(plotNSE),plotNSE,'-or');
            %plot(1:length(plotNLE),plotNLE,'-xg');
            %plot(1:length(plotNSE),plotNSE,'*');
        %end
    end
    legend(kcc2pa_fields_num);
    set(gca,'XTick',linspace(1,length(inhib_fields),length(inhib_fields)));
    set(gca,'XTickLabel',inhib_fields_num);
    xlabel('inhibition strength (x baseline conductance)');
    ylabel('"Normalisation" index');
    result = Nindex;    
end
function result = Data_2016_05_18_Diam()
    %not done
    whos
    
   %% diam comparison
    foldername = 'Data-2016-05-18(Diam)';
    keeps = {'proximal_[^K][^C]+[^2][^_].*'; 'internrn_distal_[^K][^C]+[^2][^_].*';'proximal_KCC2*';'internrn_distal_KCC2*'};
    for k = 1:length(keeps)
        keep='proximal_[^K][^C]+[^2][^_].*';   % proximal not KCC2
        keep='proximal_KCC2*';                 % proximal with KCC2
        keep = keeps{k};
        split='InPt\d+(';                       %split on inhibition strength
        degree='_[.\d]+.*InPt';                 %degrees of diam
        old = cd('Andy+Chris');
        output = getCompare(foldername,keep,splitReg,degree,'trimBeginSplit',4,'trimEndSplit',1,'trimEndDegree',4);
        cd(old);
        [inhib_fields, inhib_fields_num] = reorderByNum(fields(output));

        for inhf = 1:length(inhib_fields)
           inhib_result = output.(inhib_fields{inhf});
           cutEnd=0;    
           [diam_fields,diam_fields_num] = reorderByNum(fields(inhib_result),cutEnd);
           for diam = 1:length(diam_fields)
               %rows: same inhibition, different diam (rf,:)
               %columns: same KCC2Pa, different inhibition (:,df)
               excitation{inhf,diam} = inhib_result.(diam_fields{diam}).Excitation;
               spikes{inhf,diam} = inhib_result.(diam_fields{diam}).Spikes;
               somaEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).Soma_EGABA;
               ldendEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).ldend_EGABA;
               bdendEGABA{inhf,diam} = inhib_result.(diam_fields{diam}).bdend_EGABA_;
           end
        end
        %%TODO: put fixed cl in matrices for nindex calc
         for inhf = 1:length(inhib_fields)
             Nindex2{inhf} = InnerNormalizationIndex(excitation{inhf},spikes(inhf,:),strcat(regexprep(keep,'_','-'),inhib_fields{inhf}));
    %         %Nindex_soma_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},somaEGABA(inhf,:),inhib_fields{inhf});
    %         %Nindex_ldend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},ldendEGABA(inhf,:),inhib_fields{inhf});
    %         %Nindex_bdend_EGBABA{inhf} = InnerNormalizationIndex(excitation{inhf},bdendEGABA(inhf,:),inhib_fields{inhf});
         end
        for diam = 1:length(diam_fields)
            %Nindex between inhibition strengths for each diam
            Nindex{diam} = InnerNormalizationIndex(excitation{diam},spikes(:,diam),strcat(regexprep(keep,'_','-'),diam_fields{diam}));
        end
        figure;
        hold on;
        for n = 1:length(Nindex)
            plotN=Nindex{n}(~isnan(Nindex{n}));
            %plotNSE=Nindex_soma_EGBABA{n}(~isnan(Nindex_soma_EGBABA{n}));
            %plotNLE=Nindex_ldend_EGBABA{n}(~isnan(Nindex_ldend_EGBABA{n}));
            %plotNBE=Nindex_bdend_EGBABA{n}(~isnan(Nindex_bdend_EGBABA{n}));
            %plotNE = bdendEGABA(inhf,:);
            %[plotN,plotNE] = matchNindices(plotN,plotNE);
            %if(length(plotN)==length(inhib_fields))
                plot(1:length(plotN),plotN,'-+');
                %legplot = plot(plotNE,plotN,'-+');
                %plot(1:length(plotNSE),plotNSE,'-or');
                %plot(1:length(plotNLE),plotNLE,'-xg');
                %plot(1:length(plotNSE),plotNSE,'*');
            %end
        end
        legend(diam_fields_num);
        set(gca,'XTick',linspace(1,length(inhib_fields),length(inhib_fields)));
        set(gca,'XTickLabel',inhib_fields_num);
        xlabel('inhibition strength (x baseline conductance)');
        ylabel('"Normalisation" index');
        title(regexprep(keep,'_','-'));
        result = Nindex;    
    end
end
