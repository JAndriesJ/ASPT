function M = makeTGPT(C2Dom, lambda, ord)
    M = GPT.TGPT;
    M.lambda = lambda;
    M.order = ord;
    M = M.compTGPT(C2Dom);
    M = M.getSVDtgptMat;
end