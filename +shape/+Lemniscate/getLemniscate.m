function L = getLemniscate(foci,tStep, Sensitivity)
    if nargin < 1
        fociList = shape.Lemniscate.getListOfFoci();
        foci = fociList{randi(10,1)};
        tStep = 0.005;
        Sensitivity =tStep*5; 
    elseif nargin < 2
        if numel(foci) == 1
            fociExNum = foci;
            fociList = shape.Lemniscate.getListOfFoci();
            foci = fociList{fociExNum};
        end
        tStep = 0.001;
        Sensitivity =tStep*100; 
    end
    
    L = shape.Lemniscate.Lemniscate;
    L.foci = foci;
    L = L.getEquForm;
    L = L.getEquMotion;
    L = L.getBoundaryPoints(tStep, Sensitivity);
    L = L.getC2Information;
    L = L.interpolateC2Info;
    L = L.getC2Object;
end