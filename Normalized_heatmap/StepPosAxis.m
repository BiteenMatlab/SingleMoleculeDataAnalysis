function StepPosAxis=StepPosAxis(mainfold,maskfold,varargin)
params.From_SMAUG = false;
params.wanted_state=1;
paramsnames=fieldnames(params);
if nargin>2
    for ii=1:2:nargin-2
        whichField = strcmp(paramsnames,varargin{ii});
        try
            eval(['params.' paramsnames{whichField} ' = varargin{ii+1};'])
        catch
            error([varargin{ii}, '  is not an input parameter. Check the spelling.'])
        end
    end
end

if params.From_SMAUG
    tr=getselectedtracks(mainfold,'SMAUG_stepstates', true);
    %     for i=1:length(tr)
    %         temptr=tr{i};
    %         if isempty(temptr)
    %             continue
    %         end
    %         temptr(:,2:3)=temptr(:,2:3)/49;
    %         tr{i}=temptr(temptr(:,6)==params.wanted_state,:);
    %         clear temptr
    %     end
else
    tr=getselectedtracks(mainfold,'MinLength',4,'MaxGap',2);
end
disp('Select the phasemasks');
[maskfile,maskpath]=uigetfile([maskfold,filesep,'*Mask.mat'],'multiselect','on');
TempStepPosAxis=cell(length(tr),1);
parfor i=1:length(tr)
    tempsteps=tr{i};
    if isempty(tempsteps)
        continue
    end
    maskfilename=[maskpath maskfile{i}];
    a=load(maskfilename,'PhaseMask');
    PhaseMask=a.PhaseMask;
    
    if params.From_SMAUG
        TempStepPosAxis{i}=[DistrAxis2(tempsteps,PhaseMask),tempsteps(:,6)];
        
    else
        TempStepPosAxis{i}=DistrAxis2(tempsteps,PhaseMask);
    end
    
end
StepPosAxis=TempStepPosAxis;
% figure;
%
% scatter(StepPosAxis(:,1),StepPosAxis(:,2),1,'r');
end


%% sightly change the DistrAxis Function
function Temp_PlaceInAxis=DistrAxis2(steps,mask)

% % --select the stepsize within some range------------------------
% 
% index1=1:size(steps)-1;
% index2=index1+1;
% StepSizeIndex=sqrt((steps(index2,3)-steps(index1,3)).^2+(steps(index2,2)-steps(index1,2)).^2)<2;
% steps=steps(StepSizeIndex,:);
% %--------------------------------------------------------------

% xlist=steps(:,3);
% ylist=steps(:,2);
% 
% % labellist=steps(:,10);
% roilist=steps(:,5);
xlist = [];
ylist = [];
roilist = [];
StepSzlist = [];
TrID = unique(steps(:,4));
for i=1:length(TrID)
    temptr=steps(steps(:,4)==TrID(i),:);
    fixedTrack = nan(max(temptr(:,1)),size(temptr,2));
    cur_ROI = temptr(1,5);
    fixedTrack(temptr(:,1),:) = temptr;
    fixedTrack(1:find(all(isnan(fixedTrack),2)==0,1,'first')-1,:)=[];
    tempstepSz=sqrt(sum((fixedTrack(2:end,[3,2])-fixedTrack(1:end-1,[3,2])).^2,2));
    tempstep_pos=(fixedTrack(2:end,[3,2])+fixedTrack(1:end-1,[3,2]))/2;
    gapsID=(~isnan(tempstep_pos(:,1)))&(~isnan(tempstep_pos(:,2)));
    
    xlist = [xlist; tempstep_pos(gapsID,1)];
    ylist = [ylist; tempstep_pos(gapsID,2)];
    StepSzlist = [StepSzlist; tempstepSz(gapsID)];
    roilist =[roilist; ones(sum(gapsID),1)*cur_ROI];
end

ROINum=max(max(mask));
T = regionprops('table',mask,'PixelList');
T = feretProperties(T);

Temp_PlaceInAxis=cell(ROINum,1);
for ii=1:ROINum
    
    if T.MaxFeretDiameter(ii)>120||T.MaxFeretDiameter(ii)<30
        continue
    end
    
    tempEndPoints=T.MaxFeretDiameterEndpoints{ii};
    tempEndPoints2=T.MinFeretDiameterTrianglePoints{ii};
    p1=tempEndPoints(1,:);
    p2=tempEndPoints(2,:);
    shortp1=tempEndPoints2(1,:);% need another calculation of the foot of the triangu as the projected line(height line)
    shortp2=tempEndPoints2(2,:);
    shortp3=tempEndPoints2(3,:);
    temp_xlist=xlist(roilist==ii);
    temp_ylist=ylist(roilist==ii);
    temp_StepSzlist = StepSzlist(roilist==ii);
    %     temp_labellist=labellist(roilist==ii);
    fitsnum=length(temp_xlist);
    ROIPlaceInAxis=zeros(fitsnum,2);
    for j=1:fitsnum
        q=[temp_xlist(j) temp_ylist(j)];
        p=findprojection(p1,p2,q);
        shortp4=findprojection(shortp2,shortp1,shortp3);% find the foot one the base line of the triangu
        shortp=findprojection(shortp3,shortp4,q);
        ROIPlaceInAxis(j,1)=abs(p(1)-p1(1))/abs(p2(1)-p1(1));
        if shortp4(1)~=shortp3(1)
            ROIPlaceInAxis(j,2)=abs(shortp(1)-shortp3(1))/abs(shortp4(1)-shortp3(1));
        else
            ROIPlaceInAxis(j,2)=abs(shortp(2)-shortp3(2))/abs(shortp4(2)-shortp3(2));
        end
    end
    %     ROIPlaceInAxis(:,3)=temp_labellist;
    Temp_PlaceInAxis{ii}=[ROIPlaceInAxis, temp_StepSzlist];
end
Temp_PlaceInAxis = cat(1, Temp_PlaceInAxis{:});
end