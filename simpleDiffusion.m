function tr = simpleDiffusion(diffusionCoefficient,integrationTime,numberOfFrames,micronsPerPixel)
% diffusionCoefficient in microns per second^2
% integration time in seconds

tr = cumsum(cat(1,[0,0],sqrt(2*diffusionCoefficient*integrationTime)*...
    randn(numberOfFrames-1,2)),1);

tr(:,4:5) = tr(:,1:2)/micronsPerPixel;
tr(:,1) = 1;
tr(:,2) = 1:numberOfFrames;
tr(:,3) = nan;
end