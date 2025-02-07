function savingVideo(frames,sims,videoLength)
    rates = round(length(frames)/videoLength);
    if rates < 1
        framerates = 1;
    else 
        framerates = rates;
    end
    v = VideoWriter(strcat(sims.pathVideos,'/',sims.objectName,'_',sims.objectType,'Video'),'MPEG-4');
    v.FrameRate = framerates;
    v.Quality = 100;
    open(v);
    writeVideo(v,frames);
    close(v);
end
