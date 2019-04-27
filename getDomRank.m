function domRank = getDomRank(dir)
    [~, nbFiles] = Utility.getFileNamesInDir(dir,'/*.jpg');
    ImageDom = shape.Imgshape('Domain.jpg', 500, 1, 1, 10);
    TGPTori =  GPT.makeTGPT(ImageDom, 0.7, 1);
    domRank = nan(1, nbFiles - 1);
    
    for iFile = 2:nbFiles

        ImageRecDom = shape.Imgshape(['DomainCandidate', num2str(iFile-1), '.jpg'], 500, 1, 1, 10);
        TGPTrec =  GPT.makeTGPT(ImageRecDom, 0.7, 1);

        %%
        domRank(iFile-1) = norm(TGPTrec.TGPTmatrix -  TGPTori.TGPTmatrix)/norm(TGPTori.TGPTmatrix);
    end
end