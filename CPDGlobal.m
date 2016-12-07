function out = CPDGlobal(tr)
%
% NAME:
%       CPDGlobal
% PURPOSE:
%       Given trajectories of diffusing molecules fit all Cumulative
%       probability distributions (CPDs) of step sizes given a particular model of
%       diffusion all at once. This is a multi-domain fit, where several of
%       the  domains are the squared step size domains of the CPDs, and the
%       other is the time lag domain typically used in MSD fits to the
%       model of diffusion. 
% CATEGORY:
%       Data fitting
% CALLING SEQUENCE:
%       outStruct = CPDGlobal(tr);
% INPUTS:
%       tr:             is a cell array with each element being n X 4,
%                       where n is the number of frames. The first column
%                       is the trajectory id number (integers), the second
%                       column is the time step id (integers), and the 3rd
%                       and 4th columns are the x and y positions of the
%                       trajectory
%       
% OUTPUTS:
%       out:            Matlab structure used inside CPDGlobal to construct
%                       the required fitting function.
% USAGE:
%       Configure the analysis parameters below to correspond to your
%       particular experimental setup.
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 11/16.
% NOTES:
%       This code 'CPDGlobal.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.

%% Default analysis parameters
anProp.nMobile = 1;     % number of diffusive populations
anProp.tFrame = .0104;  % camera integration time in seconds
anProp.pixSize = .049;  % pixel size in microns
anProp.minTau = 1;      % minimum time lag in frames
anProp.maxTau = 5;      % maximum time lag in frames
anProp.overBool = 0;    % use overlapping or non-overlapping displacements?
anProp.bootNum = 300;   % number of bootstraps

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% linearized cell arrays of similar dimension
linCell=@(x)cat(1,x{:});

% disable fitting routine textual output
opts = optimset('Display','off');

%% calculate and collect all the squared step sizes for each time lag considered
for kk=1:numel(tr)
    if isempty(tr{kk})
        continue
    end
    trackNums = unique(tr{kk}(:,1))';
    
    if ~isempty(trackNums)
        for ii = trackNums
            tracks = tr{kk}(tr{kk}(:,1)==ii,[2,3,4]);
            
            % fill in the time holes with nans
            fixedTrack = nan(max(tracks(:,1)),size(tracks,2));
            fixedTrack(tracks(:,1),:) = tracks;
            
            % remove leading nans
            fixedTrack(1:find(all(isnan(fixedTrack),2)==0,1,'first')-1,:) = [];
            
            nLocs = size(fixedTrack,1);
            for jj=1:anProp.maxTau
                if anProp.overBool      % overlapping displacements
                    indvec1=jj+1:nLocs;
                    indvec2=1:nLocs-jj;
                elseif ~anProp.overBool % non-overlapping displacements
                    indvec2=1:jj:nLocs;
                    indvec1=indvec2(1:end-1);
                    indvec2=indvec2(2:end);
                end
                
                % calculate squared step sizes
                allSqSteps{kk,ii,jj}=nansum( (fixedTrack(indvec1,[2,3]) - ...
                    fixedTrack(indvec2,[2,3])).^2, 2);
            end
        end
    end
end

%% compile the cumulative probability distributions for each time lag
sqSteps=cell(anProp.maxTau,1);
for ii=1:anProp.maxTau
    wSteps = cat(1,allSqSteps{:,:,ii});
    sqSteps{ii}=sort(wSteps(wSteps > eps));     % nansum puts zeros where there were nans
end
oRanks=cellfun(@(x)linspace(0,1,numel(x))',sqSteps,'uniformoutput',0);
sqSteps = cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
nSteps = cellfun(@numel,sqSteps,'uniformoutput',0);

%% fitting function selection
funFinds = cpdFunFinder(anProp.nMobile);
cpdFun = funFinds.cpdFun;
msdFun = funFinds.msdFun;
pStart = funFinds.pStart;
bounds = funFinds.bounds;
dID = funFinds.dID;
aID = funFinds.aID;

% uninformed amplitude guesses
pStart{1}(aID) = 1/(numel(aID)+1);

%% GLOBAL FITTING
fHandle=@(p,tau,sqSteps,ranks)linCell(...
    cellfun(@(x,y)x-y,...
    cellfun(@(x,y)cpdFun(x,y,p),...
    sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0),...
    ranks,'uniformoutput',0));
eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
    sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0);

% fit to original data
y=sqSteps;
r=oRanks;

% time lag domain in seconds
tau = (1:anProp.maxTau)'*anProp.tFrame;

% fitting to the original data
[fP_nB,~,r_nB] = lsqnonlin(@(p)fHandle(...
    p,tau(anProp.minTau:end),y(anProp.minTau:end),r(anProp.minTau:end)),...
    pStart{1},bounds{1},bounds{2},opts);

% split up the residuals into into 1 for each cpd curve
rTemp = r_nB;
for ii = anProp.minTau:anProp.maxTau
    residCell{ii} = rTemp(1:nSteps{ii});
    rTemp(1:nSteps{ii}) = [];
end

% bootstrap
fP_B = zeros(anProp.bootNum,numel(pStart{1}));
parfor kk = 1:anProp.bootNum
    % resamples with replacement
    y1 = cellfun(@(x,y)sort(x(randsample(y,y,1))),sqSteps,nSteps,'uniformoutput',0);
    r1 = oRanks;
    
    % fitting to the bootstrapped data
    fP_B(kk,:) = lsqnonlin(@(p)fHandle(...
        p,tau(anProp.minTau:end),y1(anProp.minTau:end),r1(anProp.minTau:end)),...
        pStart{1},bounds{1},bounds{2},opts);
end

%% plot results
c = colormap('lines');
c = c(1:7,:);

if anProp.nMobile > 1
    nPlots = 3;
else
    nPlots = 2;
end

subplot(nPlots,1,1)
for ii = 1:numel(dID)
    h=histogram(fP_B(:,dID(ii)),'normalization','probability','displaystyle','stairs');
    set(h,'edgecolor',c(ii,:))
    hold on
end
title('Bootstrapped diffusion coefficients')
set(gca,'xscale','log'); hold off

subplot(nPlots,1,2)
for ii = anProp.minTau:anProp.maxTau
    plot(sqSteps{ii},residCell{ii}); hold all
end
title('Original data CPD residuals');
hold off

if anProp.nMobile > 1
    subplot(nPlots,1,3)
    lastAmp = ones(anProp.bootNum,1);
    for ii = 1:numel(aID)
        h=histogram(fP_B(:,aID(ii)),'normalization','probability','displaystyle','stairs');
        set(h,'edgecolor',c(ii,:))
        hold on
        
        lastAmp = lastAmp-fP_B(:,aID(ii));
    end
    h=histogram(lastAmp,'normalization','probability','displaystyle','stairs');
    title('Bootstrapped population amplitudes')
    hold off
end

% output files
out.fittedParameters_nB = fP_nB;
out.cpdResiduals = residCell;
out.sqSteps = sqSteps;
out.dID = dID;
out.aID = aID;
end