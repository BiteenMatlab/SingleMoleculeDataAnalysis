function outStruct = cpdFunFinder(nMobile)
%
% NAME:
%       cpdFunFinder
% PURPOSE:
%       Based on the number of diffusive terms, determine the fitting
%       function structure required for use in the CPDGlobal. code.
% CATEGORY:
%       Data fitting
% CALLING SEQUENCE:
%       outStruct = cpdFunFinder(nMobile);
% INPUTS:
%       nMobile:        Integer number of expected diffusive components
%       
% OUTPUTS:
%       outStruct:      Matlab structure used inside CPDGlobal to construct
%                       the required fitting function.

% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 11/16.
% NOTES:
%       This code 'cpdFunFinder.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.


% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% starting values
pStart = [.9,0,.1, .01, .0025, .00001, .2, .2, .2, .2];

% bounds for the fit
LB=-inf(1,numel(pStart));
LB(7:10) = 0;
UB=inf(1,numel(pStart));
UB(7:10) = 1;

switch nMobile
    case 1
        msdFun=@(tau,p) ...
            m2(tau,p([1,2]));
        cpdFun=@(x,y,p)1-...
            c2(x,y(1),1);
        pID = 1:2;
        
    case 2
        msdFun=@(tau,p)cat(2,...
            m2(tau,p([1,2])),...
            m2(tau,p([3,2])));
        cpdFun=@(x,y,p)1-...
            c2(x,y(1),p(4))-...
            c2(x,y(2),1-p(4));
        pID = [1:3,7];
        
    case 3
        msdFun=@(tau,p)cat(2,...
            m2(tau,p([1,2])),...
            m2(tau,p([3,2])),...
            m2(tau,p([4,2])));
        cpdFun=@(x,y,p)1-...
            c2(x,y(1),p(5))-...
            c2(x,y(2),p(6))-...
            c2(x,y(3),1-p(5)-p(6));
        pID = [1:4,7:8];
        
    case 4
        msdFun=@(tau,p)cat(2,...
            m2(tau,p([1,2])),...
            m2(tau,p([3,2])),...
            m2(tau,p([4,2])),...
            m2(tau,p([5,2])));
        cpdFun=@(x,y,p)1-...
            c2(x,y(1),p(6))-...
            c2(x,y(2),p(7))-...
            c2(x,y(3),p(8))-...
            c2(x,y(4),1-p(6)-p(7)-p(8));
        pID = [1:5,7:9];
        
    case 5
        msdFun=@(tau,p)cat(2,...
            m2(tau,p([1,2])),...
            m2(tau,p([3,2])),...
            m2(tau,p([4,2])),...
            m2(tau,p([5,2])),...
            m2(tau,p([6,2])));
        cpdFun=@(x,y,p)1-...
            c2(x,y(1),p(7))-...
            c2(x,y(2),p(8))-...
            c2(x,y(3),p(9))-...
            c2(x,y(4),p(10))-...
            c2(x,y(5),1-p(7)-p(8)-p(9)-p(10));
        pID = [1:10];
end

pStart={pStart(pID)};
bounds=[{LB(pID)},{UB(pID)}];
dID = find(ismember(pID,[1,3:6]));
aID = find(ismember(pID,7:10));

outStruct.cpdFun = cpdFun;
outStruct.msdFun = msdFun;
outStruct.pStart = pStart;
outStruct.bounds = bounds;
outStruct.dID = dID;
outStruct.aID = aID;
end