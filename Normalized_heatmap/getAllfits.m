function ALLGOODFITS=getAllfits(mainfold)

disp('Select the fit results');
[fitfile,fitpath]=uigetfile([mainfold,filesep,'*fits.mat'],'multiselect','on');
if ~iscell(fitfile)
    fitfilename=[fitpath fitfile];
    load(fitfilename,'fits');
    rows=fits.row(fits.goodfit);
    cols=fits.col(fits.goodfit);
    ROIs=fits.roinum(fits.goodfit);
    Goodfits=[rows cols ROIs];
    ALLGOODFITS=Goodfits;
else
    L=length(fitfile);
    ALLGOODFITS=cell(L,1);
    for i=1:L
        fitfilename=[fitpath fitfile{i}];
        try
            load(fitfilename,'fits')
            rows=fits.row(fits.goodfit);
            cols=fits.col(fits.goodfit);
            ROIs=fits.roinum(fits.goodfit);
            Goodfits=[rows cols ROIs];
        catch
            continue
        end
        if isempty(Goodfits)
            continue
        end
        ALLGOODFITS{i}=Goodfits;
        clear GoodTracks tracks
    end
end
end
