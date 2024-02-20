function savingVideo(frames,params,title,videoLength)
    rates = round(length(frames)/videoLength);
    if rates < 1
        framerates = 1;
    else 
        framerates = rates;
    end
    v = VideoWriter([params.path2datafiles,params.file_tag,'/',params.videoFolderName,'/',title],'MPEG-4');
    v.FrameRate = framerates;
    v.Quality = 100;
    open(v);
    writeVideo(v,frames);
    close(v);
end
