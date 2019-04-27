# Introduction
This project is a continuation of the paper work done in the paper [Identification of an algebraic domain in two dimensions from a finite number of its generalized polarization tensors](https://hal.archives-ouvertes.fr/hal-01827232/document). In it we show how one can from finitely many GPTs recover a polynomial which has level set describing the original domain. Spured on by this results we created the following paper [RECONSTRUCTION OF DOMAINS WITH ALGEBRAIC BOUNDARIES FROM
GENERALIZED POLARIZATION TENSORS]() which goes the extra mile and recovers the domain from the the polynomail level set.

## Installation
The code is pretty much stand alone except for the  [MatlabBGL toolbox](https://ch.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl) which is too large to upload here. This tool box in neccesary for the graph generation and circuit finding part of the algorithm.

## Walkthrough
The best way to see what the codes does is to walk through and example.
### Domain generation
```
nbPoints = 2000;
lambda = 0.7;
cD = shape.CompositeDoms.AddCompDom(1, nbPoints);
```
Here a domain is generated with the following information:
* C2 boundary information i.e. boundary: points, velocities, acceleration, unit normal vectors.
* The origin is contained in the boundary.
* The domain is bounded 
* The degree of the domain is known or easilly computed.

### Tesselated Generalized Polarization Tensor (TGPT)
```
TGPTcD = GPT.makeTGPT(cD, lambda, cD.degree);
```
The TGPT is calculated for the domain. 


### Domain recovery 
```
rD = recoverDomain(TGPTcD);
```
Domain recovery is done on the TGPT and is contained in the above wrapper which calls the following function.
```
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
```
Unpacking the function:
* The parameters *BifTol*, *RadInc* and *RadIni* are tunning parameters and are best illusterated in Figure...
* The first step is the creation of a object
* Next is to recover the coefficeints of the polynomial as the singular vector of the TGPT
* Then the polynomial is generated symbolicaly from the coefficients.
	* From the polynomial the levelset can be obtained.
* The bifurecation points are then found.
* Around bifurcation points there is searched for segmentation points to devide the level set into arcs.
* Between segmentation points the level set is traced out using the Hamiltonian of the polynomial.
* The arcs are kept and some discarded.
* The arcs are used to construct circuits.
* The curcuits define domain candidates.


### Plots
```
hold on
plot(cD,'blue','Linewidth',3) %Boundary of true domain
rD.plotLevelSet               %Level set of recovered polynomail
hold off

rD.plotDomCand(cD)            %Plot the recovered domains
```
### Choosing the best candidate from those recovered
```
%% Export the domains as .jpg of set pixel size
rD.exportDomCandidatesJPG     % Candidate domains
exportDomJPG(Dom)             % True domains
```
```
%% Read the exported domain .jpg files as curves, compute first TGPTs and compare 
% dir = '/home/user/MatLabProjects/PADRA/Figures/Recoverd Domains JPG Images/11';
% domRank = getDomRank(dir);
```




















