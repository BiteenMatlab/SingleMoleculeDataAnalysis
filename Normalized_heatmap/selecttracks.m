function GOODtracks=selecttracks(tracks,MinLength,MaxGap)
%% this function select the tracks with lenth > 4 steps and with a maxium gap 1
track_ID=unique(tracks(:,4));
Longstep=[];
for j=1:length(track_ID)
    if sum(tracks(:,4)==track_ID(j))>MinLength-1
        Longstep=[Longstep,track_ID(j)];
    end
end
GoodGapID=[];

for i=1:length(Longstep)
    temptrackframe=tracks(tracks(:,4)==Longstep(i),1);
    if max(diff(temptrackframe))<MaxGap+2
        GoodGapID=[GoodGapID, Longstep(i)];
    end
end
temp_tr=cell(1,length(GoodGapID));
for ii=1:length(GoodGapID)
%     label=TrackStepLabel(selectedtr);
    temp_tr{ii}=tracks(tracks(:,4)==GoodGapID(ii),:);
end
GOODtracks=cat(1,temp_tr{:});
end
