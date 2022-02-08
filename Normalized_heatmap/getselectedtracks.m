function tr=getselectedtracks(mainfold,varargin)

params.MinLength=4;%% minmum track length you want to select
params.MaxGap=5;%% maxium gap between two frame you could get
params.manualselect= true;
params.SMAUG_stepstates=false;
params.Masterfit=false;
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
if params.SMAUG_stepstates
    disp('Select the SMAUG results');
    [smaugfile,smaugpath]=uigetfile([mainfold,filesep,'*.mat'],'multiselect','off');
    smaugresult=[smaugpath smaugfile];
    load(smaugresult,'Steps','out','Sample');
    Stat= MeanCalculate(out);
    [~,D_rank]=sort(Stat.D);
    label=round(mean(Sample.LabelSaves,1))';
    D_sorted_label=zeros(length(label),1);
    for ii=1:length(D_rank)
        D_sorted_label(label==D_rank(ii))=ii;
    end
    Steps=[Steps D_sorted_label];
    fileID=unique(Steps(:,6));
    L=length(unique(Steps(:,6)));
    tr=cell(L,1);
% no need to select good track, use all steps from SMAUG
%     for i=1:L
%         tracks=Steps(Steps(:,6)==fileID(i),[10 8 9 7 5 11]);
%         GoodTracks=selecttracks(tracks,params.MinLength,params.MaxGap);
%         if isempty(GoodTracks)
%             continue
%         end
%         tr{i}=GoodTracks;
%         clear GoodTracks
%     end
 for i=1:L
     tracks=Steps(Steps(:,6)==fileID(i),[10 8 9 7]);% 5 12]);
     tracks(:,2:3)=tracks(:,2:3)/49;
     tr{i}=tracks;
     clear tracks
 end
else
    if params.manualselect==true
        disp('Select the fit results');
        if params.Masterfit
       [fitfile,fitpath]=uigetfile([mainfold,filesep,'*analysis.mat'],'multiselect','on');
       if ~iscell(fitfile); fitfile={fitfile}; end             
       L=length(fitfile);
            tr=cell(L,1);
            for i=1:L
                fitfilename=[fitpath fitfile{i}];   
            try 
                load(fitfilename,'trackfile')
                tracks=trackfile(:,[2 4 5 1 13]);
                GoodTracks=selecttracks(tracks,params.MinLength,params.MaxGap);
            catch
                continue
            end
                if isempty(GoodTracks)
                    continue
                end
                tr{i}=GoodTracks;
                clear GoodTracks tracks
            end
        else
        [fitfile,fitpath]=uigetfile([mainfold,filesep,'*.mat'],'multiselect','on');
        if ~iscell(fitfile)
            fitfilename=[fitpath fitfile];
            load(fitfilename,'tracks');
            GoodTracks=selecttracks(tracks,params.MinLength,params.MaxGap);
            tr=GoodTracks;
        else
            L=length(fitfile);
            tr=cell(L,1);
            for i=1:L
                
                fitfilename=[fitpath fitfile{i}];
                
            try 
                load(fitfilename,'tracks')            
                GoodTracks=selecttracks(tracks,params.MinLength,params.MaxGap);
            catch
                continue
            end
                if isempty(GoodTracks)
                    continue
                end
                tr{i}=GoodTracks;
                clear GoodTracks tracks
            end
        end
        end
    else
        foldprop=dir(strcat(mainfold,'\*','fits.mat'));
        L=length(foldprop);
        tr=cell(L,1);
        for i=1:L
            
            fitfilename=strcat(mainfold,'\',foldprop(i).name);
            try 
                load(fitfilename,'tracks')            
                GoodTracks=selecttracks(tracks,params.MinLength,params.MaxGap);
            catch
                continue
            end

            if isempty(GoodTracks)
                continue
            end
            tr{i}=GoodTracks;
            clear GoodTracks tracks
        end
    end
end
end