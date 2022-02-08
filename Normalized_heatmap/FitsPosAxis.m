function fitsPosAxis=FitsPosAxis(mainfold,maskfold)

ALLGOODFITS=getAllfits(mainfold);
disp('Select the phasemasks');
[maskfile,maskpath]=uigetfile([maskfold,filesep,'*Mask.mat'],'multiselect','on');
TempFitsPosAxis=cell(length(ALLGOODFITS),1);
AllMask=cell(length(maskfile),1);
for i=1:length(maskfile)
    tempfilename=[maskpath maskfile{i}];
    load(tempfilename, 'PhaseMask');
    AllMask{i}=PhaseMask;
end
parfor i=1:length(ALLGOODFITS)
% for i=1:length(ALLGOODFITS)
     tempfits=ALLGOODFITS{i};
    if isempty(tempfits)
        continue
    end
tempfitsposaxis=DistrAxis2(tempfits,AllMask{i});
TempFitsPosAxis{i}=[tempfitsposaxis ones(length(tempfitsposaxis),1)*i];
end
% fitPosAxis=TempFitsPosAxis(~cellfun('isempty',TempFitsPosAxis));
fitsPosAxis=TempFitsPosAxis;
% fitsPosAxis=cat(1,fitPosAxis{:});
% figure;
% 
% scatter(StepPosAxis(:,1),StepPosAxis(:,2),1,'r');
end


%% sightly change the DistrAxis Function
function Temp_PlaceInAxis=DistrAxis2(fits,mask)

xlist=fits(:,2);
ylist=fits(:,1);

% labellist=steps(:,10);
roilist=fits(:,3);


ROINum=max(max(mask));
if ROINum==0 
    return
end
T = regionprops('table',mask,'PixelList');
T = feretProperties(T);

Temp_PlaceInAxis=[];
for ii=1:ROINum
    %     if T.MaxFeretDiameter(ii)>100||T.MaxFeretDiameter(ii)<30
    %         continue
    %     else
%     if T.MaxFeretDiameter(ii)<60
%         continue
%     else
        
        tempEndPoints=T.MaxFeretDiameterEndpoints{ii};
        tempEndPoints2=T.MinFeretDiameterTrianglePoints{ii};
        p1=tempEndPoints(1,:);
        p2=tempEndPoints(2,:);
        shortp1=tempEndPoints2(1,:);% need another calculation of the foot of the triangu as the projected line(height line)
        shortp2=tempEndPoints2(2,:);
        shortp3=tempEndPoints2(3,:);
        temp_xlist=xlist(roilist==ii);
        temp_ylist=ylist(roilist==ii);
        if isempty(temp_xlist)
            continue
        end
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
        Temp_PlaceInAxis=[Temp_PlaceInAxis;ROIPlaceInAxis];
%     end
end
end