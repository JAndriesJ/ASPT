% This is a DEMO script in which onw caan test the functionality of the code
clear 
clc
%% Generate lemniscate domain (not fully implimented yet)
% L = shape.Lemniscate.getLemniscate(12);
% TGPT =  GPT.makeTGPT(L.C2Obj, 0.7, L.ord);
% rD = recoverDomain(TGPT);
%% Generate composite domain
nbPoints = 2000;
lambda = 0.7;
cD = shape.CompositeDoms.AddCompDom(11, nbPoints);
%%
TGPTcD = GPT.makeTGPT(cD, lambda, cD.degree);
%%
rD = recoverDomain(TGPTcD);
%% plots
hold on
% plot(cD,'blue','Linewidth',3) %Boundary of true domain
rD.plotLevelSet               %Level set of recovered polynomail
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
% legend('True Domain','Recovered polynomial','Origin','Bifurcation points','Segmentation points')
hold off
%%
rD.plotDomCand(cD)            %Plot the recovered domains
legend('True Domain','Recovered polynomial','Location','northwest')
%% Export the domains as .jpg of set pixel size
% rD.exportDomCandidatesJPG     % Candidate domains
% exportDomJPG(Dom)             % True domains
%% Read the exported domain .jpg files as curves, compute first TGPTs and compare  
% dir = '/home/andries/MatLabProjects/TGPT2AlgDom/Figures/Recoverd Domains JPG Images/11';
% rD.getDomRank(dir);
%% saving the domain and recovered domains as image.png files.
function [] = exportDomJPG(Dom)
        figure('Position', [1 1 1 1]*(400/1.5));
        plot(Dom,'red','Linewidth',3);
        axis square
        axis([-1 1 -1 1]*1.5)
        axis off
        drawnow
        pause(3)
        saveas(gcf,['Domain','.jpg'])
        close all
end

function rD = recoverDomain(TGPT)
    BifTol = 1e-4;
    RadInc = 0.01;
    RadIni = 0.05;
    
    rD = recDom;
    rD.coefVector = TGPT.singVecs(1,:);
    rD = rD.recPoly;
    rD = rD.recBifurPoints(BifTol);
    rD = rD.recSegPoints(RadIni, RadInc);
    rD = rD.recEdgeSet;
    rD = rD.recCircuits;
    rD = rD.recDomCandidate;
end


