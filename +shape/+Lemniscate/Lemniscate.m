classdef Lemniscate

    properties
        foci   = [];
        radiusSqr = 1;
        ord = 2;
        
        EquForm;
        EquMotion;
        
        OrbitTime
        Points
        C2Info
        %Velocity
        %Normals
        %Acceleration
        C2IntInfo
        %Velocity
        %Normals
        %Acceleration
        C2Obj

    end
    
    methods
      function obj = getEquForm(obj)
            syms  x y
            Foci = obj.foci;
            RadiusSqr = shape.Lemniscate.Lemniscate.getRadius(Foci); 
            
            nbFoci = length(Foci);
            XTempVec = (x.*ones(1, nbFoci) - Foci(1,:)).^2;
            YTempVec = (y.*ones(1, nbFoci) - Foci(2,:)).^2;

            Eqn.syms  = prod( XTempVec + YTempVec ) - RadiusSqr^(nbFoci);
            Eqn.func  = matlabFunction(Eqn.syms);
            
            obj.ord       = 2*nbFoci; 
            obj.radiusSqr = RadiusSqr;
            obj.EquForm   = Eqn;
      end
      
      function obj = getEquMotion(obj)
            SymEqu = obj.EquForm.syms;
            Hamiltonian = shape.Lemniscate.Lemniscate.getHamiltonian(SymEqu); 
            
            obj.EquMotion = Hamiltonian;
      end
      
      function obj = getBoundaryPoints(obj, tStep, FineSearch)
            ODEfun = obj.EquMotion.HamFun;
            IniDir = ODEfun(0,[0,0])/norm(ODEfun(0,[0,0])); 
            antiIniDir = -IniDir;
            
            EventsFunction = shape.Lemniscate.Lemniscate.makeEventsFunction(antiIniDir, FineSearch);
            tspan = 0:tStep:10;
            
            odeOptions = odeset('Events', EventsFunction,'RelTol',1e-12, 'AbsTol',1e-12);
            
            
            [Time, BndryPnts, EventTime, ~,~] = ode113(ODEfun, tspan, [0;0], odeOptions);
            BndryPnts = BndryPnts';
            OrbTime = EventTime(2);
            
            mask = (Time < OrbTime);
            obj.Points = BndryPnts(:,mask);
            obj.OrbitTime = OrbTime;
      end
      
      function obj = getC2Information(obj)
          BndyPnts = obj.Points;
          nbPoints = length(BndyPnts);
          HamFun = obj.EquMotion.HamFun;
          

          Velocity     = zeros(2,nbPoints);
          Normals      = zeros(2,nbPoints);
          Acceleration = zeros(2,nbPoints);
          
          for iPnt = 1:nbPoints
              Velocity(:,iPnt)   = HamFun(0, [BndyPnts(1,iPnt), BndyPnts(2,iPnt)]) ;
              Normals(:,iPnt)    = [0,1;-1,0]*Velocity(:,iPnt)/norm(Velocity(:,iPnt));
%               Acceleration(:,iPnt) =1 ;
          end

          obj.C2Info.Velocity = Velocity;
          obj.C2Info.Normals = Normals;
          obj.C2Info.Acceleration = Acceleration;
      end
      
      function obj = interpolateC2Info(obj)
          BndyPnts = obj.Points;
          nbPoints = length(BndyPnts);
          
          theta  = 1:nbPoints;
          fitx   = csapi(theta(:), BndyPnts(1,:));
          fity   = csapi(theta(:), BndyPnts(2,:));
          
          dfx = fnder(fitx,1);
          dfy = fnder(fity,1);
          tx = fnval(dfx, theta);
          ty = fnval(dfy, theta);
          tvec = [reshape(tx,1,[]); reshape(ty,1,[])];
          
          rotation = [[0 1];[-1 0]] ;
          normal = rotation*tvec ;
          
          normal = normal./repmat(sqrt(normal(1,:).^2+normal(2,:).^2),2,1);
          

          ddfitx = fnder(fitx,2);
          ddfity = fnder(fity,2);
          accx   = fnval(ddfitx, theta);
          accy   = fnval(ddfity, theta);
          avec   = [reshape(accx,1,[]); reshape(accy,1,[])]; 
          
          obj.C2IntInfo.Velocity   = tvec;
          obj.C2IntInfo.Normals = normal;
          obj.C2IntInfo.Acceleration   = avec;
          
      end
      
      function obj = getC2Object(obj)
          BndyPnts     = obj.Points;
          Velocity     = obj.C2IntInfo.Velocity;
          Normals      = obj.C2IntInfo.Normals;
          Acceleration = obj.C2IntInfo.Acceleration;  
          
        obj.C2Obj = shape.C2boundary(BndyPnts , Velocity, Normals, Acceleration, [], 'Lemniscate');
      end
      
      function [] = visualize(obj)
            fci     = obj.foci;
            BndryPnts = obj.Points;
            Velocity = obj.C2Info.Velocity;
            Normals  = obj.C2Info.Normals;
            Avec     = obj.C2Info.Acceleration;
            PlotBoxSize = max(vecnorm(BndryPnts,1)); 
            
            hold on
            scatter(0,0,40,'+k')
            scatter(fci(1,:),fci(2,:),'>b','filled')
            scatter(BndryPnts(1,:), BndryPnts(2,:),20,'m')
%             quiver(BndryPnts(1,:), BndryPnts(2,:), Velocity(1,:), Velocity(2,:),'g')
%             quiver(BndryPnts(1,:), BndryPnts(2,:), Normals(1,:), Normals(2,:),'c')
%             quiver(BndryPnts(1,:), BndryPnts(2,:), Avec(1,:), Avec(2,:),'b') 
            axis([-1,1,-1,1]*PlotBoxSize)
            axis square
       end
   
     end        
        methods (Static)   
            function RadiusSqr = getRadius(foci)
                Sqrs  = (foci.^2);
                radii = sum(Sqrs,1);
                RadiusSqr = prod(radii);
            end
            
            function Hamiltonian = getHamiltonian(SymEqu)
                syms x y t
                dPdx = diff(SymEqu,x);
                dPdy = diff(SymEqu,y);
                GradPoly   = [dPdx; dPdy];
                Hamiltonian.GradPoly = GradPoly;
                
                RotMat   = [0,-1;1,0];
                HamSyms = RotMat*GradPoly;
                
                HamFun = matlabFunction(HamSyms,'Vars',[t x y]);
                HamFun = @(t,z) HamFun(t,z(1),z(2));
                Hamiltonian.HamFun = HamFun;
                
                d2Pdx2  = diff(dPdx, x);
                d2Pdxdy = diff(dPdx, y);
                d2Pdydx = diff(dPdy, x);
                d2Pdy2  = diff(dPdy, y);
                HessSyms = [d2Pdx2,d2Pdydx;d2Pdxdy,d2Pdy2];
                HessFun  = matlabFunction(HessSyms,'Vars',[t x y]);
                Hamiltonian.HessFun = HessFun;
            end
            
            function EventsFunction = makeEventsFunction(antiIniDir, FineSearch)
                function [value, isterminal, direction] = EventsFcn(~, Points)
                    value      =   (norm(Points) < FineSearch)*(dot(Points, antiIniDir)  < 0);
                    isterminal =   0;
                    direction  =   1;
                end
                EventsFunction = @EventsFcn ;
            end
            
    end
end