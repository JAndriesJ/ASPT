classdef recDom
    properties
        %poly
        coefVector
        polynomial
        degree
        %splits
        bifurPoints
        segmePoints
        BifPntRadi
        SegBifRel
        %structure
        EdgeSet
        OriginEdge
        Circuits
        %Final result
        DomainCandidates
        DomainRank
    end
    
    methods
        function obj = recPoly(obj)
            coefVec = obj.coefVector;
            obj.degree = recDom.getDegree(coefVec);
            
            syms x y
            % acceptable singVec sizes
            aclengs = zeros(1,10);
            lenCoefVec = length(coefVec);
            xExpo = lenCoefVec;
            yExpo = lenCoefVec;
            count = 0 ;
            for iter1 = 2:10
                aclengs(iter1) = iter1 + aclengs(iter1-1);
                for iter2 = 0:(iter1-1)
                    count = count + 1;
                    xExpo(count) = iter1-iter2-1;
                    yExpo(count) = iter2;
                end
            end
            
            if ~any(lenCoefVec == aclengs)|| size(coefVec,1) ~= 1
                error('Not appropriate coefficient vector size')
            else
                xExpo = xExpo(1:lenCoefVec);
                yExpo = yExpo(1:lenCoefVec);
                RecPoly = sum((x.^xExpo).*(y.^yExpo).*(coefVec));
            end
            obj.polynomial = RecPoly;
        end
        
        function obj = recBifurPoints(obj,BifTol)
            searchGrid = 2;
            finess = 0.4;
            interBifDist =  0.3;
            
            poly = obj.polynomial;
            
            RFun = recDom.findRootFuncion(poly);
            sMesh = recDom.getSearchMesh(searchGrid, finess);
            
            options = optimset('Display', 'off', 'MaxIter', 20, 'Algorithm', 'Levenberg-Marquardt');
            BifPnts = recDom.getBifurPoints(RFun, sMesh, BifTol, options, interBifDist);
            
            obj.bifurPoints = BifPnts;
        end
        
        function obj = recSegPoints(obj, RadIni, RadInc) 
            poly = obj.polynomial;
            BifPnts = obj.bifurPoints;
            
            minDistCtrl = 1000 ;
            
            [SegPnts, Radii, SegBifRelCel] =  recDom.getSegPnts(poly, BifPnts, minDistCtrl, RadIni, RadInc);
            BifPntRadii = Radii ;
            SegPnts = [SegPnts ; 1:size(SegPnts,2)]; % Label segmentation points
            
            obj.segmePoints = SegPnts;
            obj.BifPntRadi = BifPntRadii;
            obj.SegBifRel = SegBifRelCel;
        end
        
        function obj = recEdgeSet(obj)
            poly         = obj.polynomial;
            SegPts       = obj.segmePoints;
            BifRadii     = obj.BifPntRadi;
            SegBifRelCel = obj.SegBifRel;
            
            [EdgSet, OriEdge] =  ... 
                recDom.getEdgeRel(poly, SegPts, BifRadii, SegBifRelCel);
            
            obj.EdgeSet = EdgSet;
            obj.OriginEdge = OriEdge;
        end
        
        function obj = recCircuits(obj)
           EdgSet  = obj.EdgeSet;
           OriEdge = obj.OriginEdge;
           SgBfRl  = obj.SegBifRel;
           
           addpath('matlab_bgl')
           ESet = EdgSet(:,1:2) ;
           G =  digraph(ESet(:,1), ESet(:,2)) ;
           [~,CirSet] = find_elem_circuits(G.adjacency) ;
           rmpath('matlab_bgl');
           
           CirSet = recDom.removeShortLongCircuits(CirSet, SgBfRl);
           
           CirSet = recDom.removeTrivCircuits(CirSet, EdgSet);
           
           CirSet = recDom.removeNonOriginCircuits(CirSet, OriEdge);
           
           CirSet = recDom.removeDupliCircuits(CirSet);
           
           CirSet = recDom.removePinchCircuit(CirSet, SgBfRl);
    
           obj.Circuits = CirSet;
        end
        
        function obj = recDomCandidate(obj)
            poly     = obj.polynomial;
            deg      = obj.degree;
            SegPts   = obj.segmePoints;
            BifRadii = obj.BifPntRadi;
            Edgs     = obj.EdgeSet;
            Cirs     = obj.Circuits;
            
            nbCir = length(Cirs);
            for iCir = 1:nbCir
                Cir = Cirs{iCir};
                DomCandidate = recDom.getDomCandidate(poly, SegPts, BifRadii, Edgs, Cir);
                
                % Shift the beginning to the origin
                PointNorms = vecnorm(DomCandidate');
                [~,indZero] = min(PointNorms);
                shift = length(DomCandidate) - indZero;
                DomCandidate = circshift(DomCandidate', shift, 2);
                
                theta = 1:length(DomCandidate);
                [points, tvec, avec, normal] = ...
                    shape.C2boundary.boundary_vec_interpl(DomCandidate, theta);
                
                rDom = shape.C2boundary(points, tvec, avec, normal, [0;0], 'recDom');
                
                obj.DomainCandidates{iCir} = rDom;
            end
            
            
        end
        
        function obj = rankDomCandidates(dir)
            obj.DomainRank = getDomRank(dir);
        end
        
        function [] = plotLevelSet(obj)
            hold on
            % Polynomial level-set
            poly = obj.polynomial;
            recDom.plotPoly(poly)
            scatter(0,0,150,'k','+')
            
            % Bifurcation points
            recDom.plotBif(obj.bifurPoints);

            % Segmentation points
            recDom.plotSeg(obj.segmePoints)


            axis square
            hold off
        end
        
        function [] = plotDomCand(obj, TrueDom) 
            if (nargin == 1)
                TrueDom.points = [0;0];
            end
            DomCand = obj.DomainCandidates;
            DomRank = obj.DomainRank;
            % Domain candidates
            recDom.plotDomCands(DomCand, TrueDom, DomRank)
        end
        
        function [] = exportDomCandidatesJPG(obj)
            DomCands = obj.DomainCandidates;
            nbDom = length(DomCands);
            for iDom = 1:nbDom
                figure('Position', [1 1 1 1]*(400/1.5));
                plot(DomCands{iDom},'Linewidth',4);
                axis square
                axis([-1 1 -1 1]*1.5)
                axis off
                pause(3)
                saveas(gcf,['DomainCandidate', num2str(iDom),'.jpg'])
                close all
            end
        end
        
    end
    
    methods(Static)
        %% poly
        % return degree of domain bassed on recovered coefficients
        function deg = getDegree(Coef)
            maskDeg = (length(Coef) == cumsum(2:9));
            degVec = (1:8);
            deg = degVec(maskDeg);
        end
        
        %% Recover bifurcation points
        % main search
        function BifPnts = getBifurPoints(RFun, sMesh, BifTol, options, interBifDist)
            BifPnts = [] ;
            isLessThan2 = (size(BifPnts,2) < 2);
            isPatient = 1;
            count = 1;
            while isLessThan2 && isPatient
                
                BifPnts = recDom.searchBifurPoints(RFun, sMesh, BifTol, options, interBifDist);
                BifPnts = recDom.removeOrigin(BifPnts);
                BifTol  = recDom.incSearchTol(BifTol);
                
                count = count +1;
                isLessThan2 = (size(BifPnts,2) < 2);
                isPatient = count < 10;
            end
        end
        
        function BifPnts = searchBifurPoints(RFun, sMesh, BifSearchTol, options, interBifDist)
            count = 0;
            BifPnts    = zeros(10,2) ;
            for iter = 1 : size(sMesh,2)
                searchPoint = [sMesh(1,iter),sMesh(2,iter)];
                candidateBifurPoint = fsolve(RFun , searchPoint, options);
                
                isBifur = norm(RFun(candidateBifurPoint)) <  BifSearchTol;
                isDestinct = all(vecnorm(BifPnts'-candidateBifurPoint') > interBifDist);
                
                if isBifur && isDestinct
                    count = count + 1;
                    BifPnts(count,:) = candidateBifurPoint;
                end
            end
        end
        
        function NoOrigBifPnts = removeOrigin(BifPnts)
            maskNonZero = any(BifPnts, 2);
            NoOrigBifPnts = BifPnts(maskNonZero, :)';
        end
        
        function BifSearchTol = incSearchTol(BifSearchTol)
            toloranceMultiplier  = 2;
            BifSearchTol = BifSearchTol*toloranceMultiplier ;
        end
        
        %Grid parameters
        function searchMesh = getSearchMesh(searchGrid, finess)
            searchVec = -searchGrid:finess:searchGrid;
            searchMesh = combvec(searchVec,searchVec);
        end
        % Root function
        function RootFunction = findRootFuncion(poly)
            deriv = recDom.derivativez(poly);
            RootFunctionTemp = matlabFunction(deriv) ;
            RootFunction = @(z) RootFunctionTemp(z(1),z(2));
        end
        
        function F = derivativez(poly)
            syms x y
            F(1) = poly ;
            F(2) = diff(poly,x) ;
            F(3) = diff(poly,y) ;
        end
        
        %% Recover segmentation points
        function [SegPnts, Radii, SegBifRelCel] =  getSegPnts(poly, BifurPoints, minDistCtrl, RadIni, RadInc)
            nbBifurPts    = size(BifurPoints,2);
            
            SegBifRelCel  = cell(1,nbBifurPts);
            SegPntsCELL   = cell(1,nbBifurPts);
            Radii         = zeros(1,nbBifurPts);
            
            Count = 1;
            for iBifPnt = 1:nbBifurPts
                searchCenter   = BifurPoints(1:2,iBifPnt) ;
                [ZeroAngles, Radius] = recDom.getZerosAngles(poly, searchCenter, minDistCtrl, RadIni, RadInc);
                
                SegPntsCELL{iBifPnt} =  recDom.convPolToCart(searchCenter, Radius ,ZeroAngles);
                
                Radii(iBifPnt) = Radius ;
                nbCandPts      = length(ZeroAngles);
                [SegBifRelCel{iBifPnt}, Count] = recDom.getSegBifurRelCel(nbCandPts, Count);
            end
            SegPnts = cell2mat(SegPntsCELL);
        end
        
        %%Search for zeros in a circle around a bifurcation point.
        function [ZeroAngles, Radius] = getZerosAngles(poly, center, minDistCtrl, RadIni, RadInc)
            ZeroAngles = [] ;
            Radius = RadIni ;
            while length(ZeroAngles) < 4
                Radius = Radius +RadInc;
                F = recDom.getObjFunct(poly, center, Radius) ;
                ZeroAngles =  recDom.findAllZeros(F, 0, 2*pi, minDistCtrl) ;
            end
        end
        
        function F = getObjFunct(poly, center, radius)
            syms x y
            TempFunc = matlabFunction(poly) ;
            F = @(theta) TempFunc(radius*cos(theta) + center(1,:), radius*sin(theta) + center(2,:)) ;
        end
        
        function z= findAllZeros(f, xmin, xmax, N)
            % Inputs :
            % f : function of one variable
            % [xmin - xmax] : range where f is continuous containing zeros
            % N : control of the minimum distance (xmax-xmin)/N between two zeros
            if (nargin < 4)
                N = 100;
            end
            dx = (xmax-xmin)/N;
            x2 = xmin;
            y2 = f(x2);
            z = [];
            options = optimset('Display', 'off');
            for iter = 1 : N
                x1 = x2;
                y1 = y2;
                x2 = xmin+iter*dx;
                y2 = f(x2);
                if (y1*y2 <= 0)                       % Rolle's theorem : one zeros (or more) present
                    z = [z,fzero(f,[x1,x2],options)]; % Linear approximation to guess the initial value in the [x1,x2] range.
                end
            end
        end
        
        % Utility : polar to Cartesian conversion
        function points = convPolToCart(center, Radius ,ZeroAngles)
            x = Radius*cos(ZeroAngles)+center(1);
            y = Radius*sin(ZeroAngles)+center(2);
            
            points = [x; y];
        end
        
        %
        function [SegBifurRelCel, Count] = getSegBifurRelCel(nbCandPts, Count)
            SegBifurRelCel = Count : (Count + nbCandPts - 1);
            Count = Count + nbCandPts;
        end
        
        %% Recover Edge relations
        %Find all the edges between segmentation points.
        function [EdgeSet, OriginEdge] = getEdgeRel(poly, SegPts, BifRadii, SegBifRelCel)
            StopRadius = max(BifRadii);
            TerRadius = 0.25*StopRadius;
            Hamiltonian = recDom.getHamiltonian(poly);
            
            TrivEdgeSet = recDom.getTrivEdges(SegBifRelCel);
            
            nbSegPts = size(SegPts,2);
            EdgeCell = cell(nbSegPts, 1);
            OriginEdge = [] ;
            for iSeg = 1:nbSegPts
                CurSeg = SegPts(1:2,iSeg);
                IndCurSeg = SegPts(3,iSeg);
                maskNotCurSeg = ~(1:nbSegPts == iSeg);
                OtherSigPnts = SegPts(1:2,maskNotCurSeg);
                
                
                EdgePoints    = ...
                    recDom.traceEdge(Hamiltonian, OtherSigPnts, CurSeg,  TerRadius);
                
                indStopSeg = recDom.getStopSeg(EdgePoints, SegPts(1:2,:), StopRadius);
                
                
                if isempty(OriginEdge)
                    [OriginEdge, ~] = recDom.checkIfOriEdge(EdgePoints, IndCurSeg,  indStopSeg);
                end 
                
                [isTriv]   = recDom.checkTrivEdge(IndCurSeg, indStopSeg, TrivEdgeSet);
                if  isTriv 
                    Edges = [];
                else
                    Edges = [IndCurSeg, indStopSeg, 1 ;...
                             indStopSeg, IndCurSeg, 2  ];
                end
                EdgeCell{iSeg} = Edges;
            end
            
            EdgeSet = cell2mat(EdgeCell);
            EdgeSet = [TrivEdgeSet;EdgeSet];
        end
        
        function TrivEdgeSet = getTrivEdges(SegBifRelCel)
            nbBifPnts = length(SegBifRelCel);
            TrivEdgeSet = [];
            for iBifPnt = 1:nbBifPnts
                SegGroup = combvec(SegBifRelCel{iBifPnt}, SegBifRelCel{iBifPnt});
                indSegGroup = (SegGroup(1,:) ~= SegGroup(2,:));
                TrivEdgeSet = [TrivEdgeSet, SegGroup(:, indSegGroup)];
            end
            nbEdges = length(TrivEdgeSet);
            TrivEdgeSet = [TrivEdgeSet', zeros(nbEdges,1)];
        end
        %Tracing out the polynomial between two segmentation points.
        function EdgePoints = traceEdge(Hamiltonian, OtherSigPnts, CurrentSeg,  TerRadius)

            tspan = 0:0.001:10;
            ode113options = recDom.getODEoptions(OtherSigPnts, TerRadius);
            [~,EdgePoints] = ode113(Hamiltonian, tspan, CurrentSeg, ode113options) ;
        end
        
        function odeOptions = getODEoptions(OtherSigPnts, TerRadius)
            absoluteTOL  = 1e-13;
            relaltiveTOL = 1e-13;
            maxNbSteps   = 1e-3;
            
            EventFunc = @(T, Y) recDom.EventsFunction(T, Y, OtherSigPnts, TerRadius);
            odeOptions = odeset('Events', EventFunc, 'AbsTol', absoluteTOL, 'Reltol', relaltiveTOL, 'MaxStep', maxNbSteps);
        end
        
        function [value, isterminal, direction] = EventsFunction(~, Y, OtherSigPnts ,TerRadius)
            OuterRadius = 1.5;
            
            OuterDom  = double(norm(Y) < OuterRadius);
            
            InnerDoms = (vecnorm([Y(1);Y(2)]-OtherSigPnts) < TerRadius);
            
            value = any(InnerDoms)*OuterDom; 
            isterminal = 1;
            direction = 0;
        end
         %Finding the index of the segmentation point where the trace stops
        function indStopSeg = getStopSeg(EdgePoints, TerSegments, StopRadius)
            edgeEndPoint = EdgePoints(end,:)';
            StopDists = vecnorm(edgeEndPoint - TerSegments);
            [MinDist, indStopSeg] = min(StopDists);
            isOutOfBounds = (MinDist > 2*StopRadius);
            if isOutOfBounds
                indStopSeg = [] ;
            end
        end
        %Determine if an edge passes through the origin 
        function [OriginEdge, isTheOriginEdge] = checkIfOriEdge(Y, CurrentSegInd, StopSegInd)
            OriginEdge = [];
            isTheOriginEdge =  (min(vecnorm(Y')) < 1e-2);
            if isTheOriginEdge
                OriginEdge = [CurrentSegInd, StopSegInd, 1 ;
                              StopSegInd, CurrentSegInd, 2];
            end
        end
        % detrmine if the edge is trivial
        function [isTriv] = checkTrivEdge(IndCurSeg, indStopSeg, TrivEdgeSet)
             if isempty(indStopSeg) || (IndCurSeg == indStopSeg)
                 isTriv = true;
                 return
             end
            
             maskStart = (IndCurSeg == TrivEdgeSet(:,1));
             maskEnd = (indStopSeg == TrivEdgeSet(:,2));
             isTriv = any(maskStart & maskEnd);
        end
       
        %% Recover circuits
        %%Circuit Lengths
        % Remove circuits that are too short or too long
        function CircuitSet = removeShortLongCircuits(CircuitSet, SegBifurRelCel)
            circuitLengths = recDom.getLenCircuits(CircuitSet);
            [maskCirShort, maskCirLong] = recDom.getCircuitLenMask(circuitLengths, SegBifurRelCel);
            maskRightLength = ~(maskCirShort | maskCirLong);
            CircuitSet = CircuitSet(maskRightLength);
        end
        
        function circuitLengths = getLenCircuits(CircuitSet)
            circuitLengths = cellfun(@length,CircuitSet);
        end
        
        function [maskCirShort, maskCirLong] = getCircuitLenMask(circuitLengths, SegBifurRelCel)
            nbBifPnts = length(SegBifurRelCel);
            
            maskCirShort = (circuitLengths <= 4);
            maskCirLong = (circuitLengths >= 3*nbBifPnts);
        end
        
        %%Edge Types
        % Remove circuits whith two consecutive trivial edges
        function CircuitSet = removeTrivCircuits(CircuitSet, EdgeSet)
            EdgeTypeCirCell = recDom.getCirCellEdgeType(CircuitSet, EdgeSet);
            maskTrivCir = recDom.getMaskTrivCir(EdgeTypeCirCell);
            maskNonTrivCir = ~maskTrivCir;
            CircuitSet = CircuitSet(maskNonTrivCir);
        end
        
        function EdgeTypeCirCell = getCirCellEdgeType(CircuitSet, EdgeSet)
            nbCir = length(CircuitSet);
            TempFunct = @(C) recDom.getCirEdgeType(C, EdgeSet);
            
            EdgeTypeCirCell = cell(1,nbCir);
            for iCir = 1:nbCir
                EdgeTypeCirCell{iCir} = TempFunct(CircuitSet{iCir});
            end
        end
        
        function EdgeTypeCir = getCirEdgeType(Circuit, EdgeSet)
            cirEdges = recDom.getCirEdges(Circuit);
            [Lia, Locb] = ismember(EdgeSet(:,1:2), cirEdges, 'rows');
            EdgeTypeUnOrd = EdgeSet(Lia,3);
            OrdFound = Locb(Locb ~= 0);
            [~,Perm] = sort(OrdFound);
            EdgeTypeCir = EdgeTypeUnOrd(Perm)';
        end
        
        function cirEdges = getCirEdges(Circuit)
            cirEdges = [Circuit(1:(end-1))', Circuit(2:(end))'];
        end
        
        % masks
        function maskTrivCir = getMaskTrivCir(EdgeTypeCirCell)
            maskTrivCir = cellfun(@recDom.getIsTrivCir, EdgeTypeCirCell);
        end
        
        function isTrviCir = getIsTrivCir(EdgeTypeCir)
            typePairs = [EdgeTypeCir(1:end);  circshift(EdgeTypeCir,1)];
            maskZeroEdgePair = (sum(typePairs,1)==0);
            isTrviCir = any(maskZeroEdgePair);
        end
        
        %%Circuits with origin edge
        %Remove circuits that do not contain the origin edge
        function CircuitSet = removeNonOriginCircuits(CircuitSet, OriginEdge)
            maskContOriEdge = recDom.getMaskContOriEdge(OriginEdge, CircuitSet);
            CircuitSet = CircuitSet(maskContOriEdge);
        end
        
        function maskContOriEdge = getMaskContOriEdge(OriginEdge, CircuitSet)
            checkOriEdgeFun = @(C) recDom.getIsContOriEdge(OriginEdge, C);
            maskContOriEdge = cellfun(checkOriEdgeFun, CircuitSet);
        end
        
        function isContOriEdge = getIsContOriEdge(OriginEdge, Circuit)
            cirEdges = recDom.getCirEdges(Circuit);
            [maskMem, ~] = ismember(OriginEdge(:,1:2), cirEdges, 'rows');
            isContOriEdge = any(maskMem);
        end
        
        %%Duplicate circuits
        %Remove circuits that are the same upto direction
        function CircuitSet = removeDupliCircuits(CircuitSet)
            maskDuplicate = recDom.getMaskDuplicate(CircuitSet);
            maskNonDup = ~maskDuplicate;
            CircuitSet = CircuitSet(maskNonDup);
        end
        
        function maskDuplicate = getMaskDuplicate(CircuitSet)
            nbCir = length(CircuitSet);
            maskDuplicate = false(1, nbCir);
            for iCir = nbCir:-1:2
                for iCirRest = (iCir-1):-1:1
                    isDup = recDom.getIsReverse(CircuitSet{iCir}, CircuitSet{iCirRest});
                    if isDup
                        maskDuplicate(iCir) = true;
                        break
                    end
                end
            end
        end
        
        function isReverse = getIsReverse(Circuit, CircuitAlt)
            isReverse = isequal(Circuit, fliplr(CircuitAlt));
        end
        
        % pinch points
        function CircuitSet = removePinchCircuit(CircuitSet, SegBifRel)
            maskPinch = recDom.getMaskPinch(CircuitSet, SegBifRel);
            CircuitSet = CircuitSet(maskPinch);
        end
        
        function maskPinch = getMaskPinch(CircuitSet, SegBifRel)
            nbCir = length(CircuitSet);
            maskPinch = zeros(1,nbCir);
            for iCir = 1:nbCir
                Circuit = CircuitSet{iCir};
                bifCir = recDom.getBifurCircuit(Circuit, SegBifRel);
                maskPinch(iCir) = (2*length(unique(bifCir)) == length(bifCir));
            end
            maskPinch = logical(maskPinch);
        end     
        
        function bifCircuit = getBifurCircuit(Circuit, SegBifRel)
            nbSeg = length(Circuit)-1;
            nbBif = length( SegBifRel);
            indBif = 1:nbBif;
            repSeg = cell(1,nbBif);
            bifCircuit = nan(1,nbSeg);
            for iSeg = 1:nbSeg
                [repSeg{:}] = deal(Circuit(iSeg));
                maskBif = cellfun(@ismember, repSeg, SegBifRel);
                bifCircuit(iSeg) = indBif(maskBif);
            end
        end
        
        %% Recover Domains
        function DomCandidate = getDomCandidate(poly, SegPts, BifRadii, Edges, Circuit)
            StopRadius = max(BifRadii);
            TerRadius = 0.35*StopRadius;
            
            nbEdges    = length(Circuit);
            EdgePoints = cell(1, nbEdges-1);
            for iEdge = 1:(nbEdges-1)
                [CurEdge, CurSeg, NextSeg ] = ...
                    recDom.loadEndPoints(SegPts, Circuit, iEdge);
                Direction = recDom.getDirection(Edges, CurEdge);
                
                if (Direction ~= 0)
                    Hamiltonian = recDom.getHamiltonian(poly, Direction);
                    EdgePoints{iEdge} = ...
                        recDom.traceEdge(Hamiltonian, NextSeg, CurSeg, TerRadius);
                end
            end
            
            maskEmptyEdge = cellfun(@isempty, EdgePoints(1:end));
            isTrivEnd = (maskEmptyEdge(end) == 0);
            
            TempEdgePoints = recDom.cirCellShift(EdgePoints, isTrivEnd);
            
            % interpolation stage looks good but does it generalize
            TempEdgePoints = recDom.interpolateTrivEdges(TempEdgePoints, nbEdges);

            %Undo the circular shift
            EdgePoints = recDom.UndoCirCellShift(TempEdgePoints, isTrivEnd);

            DomCandidate = cat(1,EdgePoints{:});
        end
        
        %Returns various points of interest
        function [CurEdge, CurSeg, NextSeg ] =...
                loadEndPoints(SegPts, Circuit, iEdge)
            
            indCurSeg = Circuit(iEdge);
            indNextSeg = Circuit(iEdge+1);
            
            CurEdge = [indCurSeg, indNextSeg];
            
            CurSeg = SegPts(1:2, indCurSeg);
            NextSeg = SegPts(1:2, indNextSeg);
        end
        
        %Gets the direction of motion for the point
        function Direction = getDirection(Edges, CurEdge)
            [Lia, ~] = ismember(Edges(:, 1:2), CurEdge, 'rows');
            Direction = Edges(Lia, 3);
        end
        
        %Interpolate all the trivial edges in a circuit
        function TempEdgePoints = interpolateTrivEdges(TempEdgePoints, nbEdges)
            for iEdge = 2:2:(nbEdges-1)
                island1 = TempEdgePoints{iEdge-1}((end-9):end,:);
                island2 = TempEdgePoints{iEdge+1}(1:10,:);
                points0 = [island1', island2'];
                
                interDist1 = mean(vecnorm(diff(island1)'));
                interDist2 = mean(vecnorm(diff(island2)'));
                
                GapDist =  norm(island1(end,:) - island2(1,:));
                theta0 = [0:interDist1:(interDist1*9), ((interDist1*9)+GapDist):interDist2:(interDist2*18 +GapDist) ];
                theta  = (interDist1*10):interDist1:((interDist1*8)+GapDist);
                Pointa = recDom.interpolatePoints(points0, theta0, theta);
                TempEdgePoints{iEdge} = Pointa';
            end
        end
        
        %Interpolat the zero direction edges.
        function interpolPoints = interpolatePoints(points0, theta0, theta)
            fx = csapi(theta0(:), points0(1,:));
            px = fnval(fx, theta);
            
            fy = csapi(theta0(:), points0(2,:));
            py = fnval(fy, theta);
            
            interpolPoints = [reshape(px,1,[]); reshape(py,1,[])];
        end
        
        %Shift the edge points in order to have non trivial edges on the..
        function TempEdgePoints = cirCellShift(EdgePoints, isTrivEnd)
            if isTrivEnd
                TempEdgePoints = [EdgePoints{end}; EdgePoints(:)]';
            else
                TempEdgePoints = [EdgePoints{:}; EdgePoints(1)]';
            end
        end
        
        
        function EdgePoints = UndoCirCellShift(TempEdgePoints, isTrivEnd)
            if isTrivEnd
                EdgePoints = TempEdgePoints(2:end);
            else
                EdgePoints = TempEdgePoints(1:end-1);
            end
        end
        
        %% Recover TGPTs
        
        %% Plot
        %Plot the polynomial levelset
        function [] = plotPoly(poly,rang)
            if nargin < 2
                rang = 1.5;
            end
            f = matlabFunction(poly) ;
            fc = fcontour( f,'r') ;
            fc.LineWidth = 2 ;
            fc.LevelList = 0 ;
            fc.MeshDensity = 500 ;
            fc.XRange = [-1,1]*rang ;
            fc.YRange = [-1,1]*rang ;
        end
        
        %Plot bifurcation points
        function [] = plotBif(bifPoints)
            noBif = isempty(bifPoints);
            if ~noBif
                bifPnts = bifPoints;
                scatter(bifPnts(1,:), bifPnts(2,:), 40, 'green','<','filled')
            end
        end
        
        %Plot segmentation points
        function [] = plotSeg(segPoints)
            noSeg = isempty(segPoints);
            if ~noSeg
                scatter(segPoints(1,:),segPoints(2,:), 40, 'magenta','o','filled')
                
                nbSeg = length(segPoints(3,:));
                SegText = cell(nbSeg,1);
                for iter = 1:nbSeg
                    SegText{iter} = num2str(segPoints(3,iter)) ;
                end
                text(segPoints(1,:),segPoints(2,:),SegText)
            end
        end
        
        %Plot domain candidates
        function [] = plotDomCands(DomCand, TrueDom, DomRank)
            if (nargin == 1)
               TrueDom.points = [0;0];
            end
            
            if ~isempty(DomCand)
                figure
                nbCandidates = length(DomCand);
                nbSubPlot = ceil(sqrt(nbCandidates));
                for icDom = 1:1%nbCandidates
                    Utility.subtightplot(1,1,icDom)
                    hold on
                    if ~isempty(DomRank)
                        text(-0.7,-0.7,['RelErr: ',num2str(round(DomRank(icDom),3))])
                    end
                    s = scatter(TrueDom.points(1,:), TrueDom.points(2,:), 15, 'blue','filled');
%                     s.MarkerFaceAlpha = 0.1;
                    plot(DomCand{icDom},'color','red','Linewidth',2);
                    hold off
                    axis([-1 1 -1 1]*1.5)
                    set(gca,'YTickLabel',[])
                    set(gca,'XTickLabel',[])
                    axis square
                end
            end
        end
        
        %% Utility
        % Used to trace out the polynomial as an equation of motion.
        function Hamiltonian = getHamiltonian(poly, direction)
            if (nargin == 1)
               dir = 1;
            else
                if (direction == 2)
                    dir = -1;
                else
                    dir = 1;
                end
            end
            
            syms x y t
            dfdx = diff(poly,x);
            dfdy = diff(poly,y);
            gradf = [dfdx; dfdy];
            
            GradPoly   = gradf/norm(gradf) ;
            RotGradPos = (dir*[0,-1;1,0])*GradPoly;
            fRotGrad= matlabFunction(RotGradPos,'Vars',[t x y]);
            Hamiltonian = @(t,z) fRotGrad(t,z(1),z(2));
        end
       
    end

end