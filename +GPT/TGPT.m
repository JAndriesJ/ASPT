classdef TGPT
    properties
        lambda
        order
        NPOpeMat
        TGPTmatrix
        singVecs
        singVals
    end  
    
    methods  
        function obj = compTGPT(obj, C2Dom)
            lamb = obj.lambda;
            ord = obj.order;
            
            npOpeMat = GPT.TGPT.getNPOpeMat(C2Dom,lamb) ;
            obj.NPOpeMat = npOpeMat ;
            
            TGPTmat = GPT.TGPT.getTGPTmat(C2Dom,npOpeMat,ord) ;
            obj.TGPTmatrix = TGPTmat ;
            
%         obj.NoisyTGPTmat = TGPTmat + Noise.*(2*rand(size(TGPTmat,1),size(TGPTmat,2))-ones(size(TGPTmat,1),size(TGPTmat,2))).*TGPTmat ;
        end
        
        function obj = getSVDtgptMat(obj)
            tgpt = obj.TGPTmatrix ;    
            [~,DiagMat,UniMat]  = svds(tgpt,1,'smallest','Tolerance',1e-10,'MaxIterations',500) ;
            sngVals = diag(DiagMat) ; 
            sngVec  = UniMat';
            
            obj.singVals = sngVals;
            obj.singVecs  = sngVec;
        end
        
    end
    methods (Static)
        
        % Generate full TGPT matrix 
        function TGPTmat = getTGPTmat(C2Dom,NPOpeMat,ord)
        indexVec = [0,cumsum(2:(2*ord+1))] ;
        TGPTmat = zeros(((2*ord+1)*(2*ord+2))/2-1,((ord+1)*(ord+2))/2-1) ;
        
        for itr1 = 1:ord
            for itr2 = itr1:ord
                if itr1 == itr2
                    % top part
                    TGPTmat((indexVec(itr1)+1):indexVec(itr1+1),(indexVec(itr2)+1):indexVec(itr2+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr2,itr1) ;
                    % bottom part
                    TGPTmat(((indexVec(itr1+ord)+1):indexVec(itr1+ord+1)),(indexVec(itr2)+1):indexVec(itr2+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr2,itr1+ord) ;
                else
                    % top part
                    TGPTmat((indexVec(itr1)+1):indexVec(itr1+1),(indexVec(itr2)+1):indexVec(itr2+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr2,itr1) ;
                    TGPTmat((indexVec(itr2)+1):indexVec(itr2+1),(indexVec(itr1)+1):indexVec(itr1+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr1,itr2) ;
                    % bottom part
                    TGPTmat(((indexVec(itr1+ord)+1):indexVec(itr1+ord+1)),(indexVec(itr2)+1):indexVec(itr2+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr2,itr1+ord) ;
                    TGPTmat(((indexVec(itr2+ord)+1):indexVec(itr2+ord+1)),(indexVec(itr1)+1):indexVec(itr1+1)) =...
                        GPT.TGPT.PTblock(C2Dom,NPOpeMat,itr1,itr2+ord) ;                    
                end               
            end
        end
        end
        
        % Generate a tessera of GPTs 
        function Mnm = PTblock(C2Dom,NPOpeMat,m,n)           
            if sum(m) == 0 || sum(n) == 0
                error('alpha and beta must be larger than zero in either coordinates')
            end

            points = C2Dom.points ;
            normal = C2Dom.normal ;
            sigma =  C2Dom.sigma  ;
            
            fun1 = @(a,b) a.*b.^(a-1);
            fun2 = @(a,b) a.^(b);

            X1  = bsxfun(fun1,(n:-1:0),points(1,:)') ;
            X2  = bsxfun(fun1,(0:n),points(2,:)') ;

            X1(~isfinite(X1))  = 0 ;
            X2(~isfinite(X2))  = 0 ;

            normalxX1 =  bsxfun(@times,X1,normal(1,:)') ;
            normalxX2 =  bsxfun(@times,X2,normal(2,:)') ;
            phi = NPOpeMat\(normalxX1+normalxX2) ;

            Y = (bsxfun(fun2,points(1,:),(m:-1:0)')).*(bsxfun(fun2,points(2,:),(0:m)')) ;
            Yxsigma = bsxfun(@times,Y,sigma) ;

            Mnm = conj(Yxsigma*phi)' ;
        end
        
        %Genrate NP-operator kernal matrix minuse the lambda
        function NPOpeMat = getNPOpeMat(C2Dom, lambda)
            nbPoints = C2Dom.nbPoints;
            NPOpeMat = lambda*eye(nbPoints) - GPT.TGPT.getNPOpeKerMat(C2Dom); 
        end

        %Genrate NP-operator kernal matrix
        function NPOpeKerMat = getNPOpeKerMat(C2Dom)
            nbPoints = C2Dom.nbPoints;
            NPOpeKerMat = zeros(nbPoints, nbPoints);
            points    = C2Dom.points;
            normal    = C2Dom.normal;
            sigma     = C2Dom.sigma;
            avec      = C2Dom.avec;
            tvec_norm = C2Dom.tvec_norm;
            
            for j = 1:nbPoints
                
                XdotY = (points(1,j)-points(1,:)).*normal(1,j)...
                    +(points(2,j) - points(2,:)).*normal(2,j);
                
                norm_XY_square = (points(1,j)-points(1,:)).^2 ...
                    +(points(2,j)-points(2,:)).^2;
                
                %Off diagonal row terms
                NPOpeKerMat(j, 1:j-1) = 1/(2*pi)* ...
                    XdotY(1:j-1)./norm_XY_square(1:j-1).*C2Dom.sigma(1:j-1);
                
                NPOpeKerMat(j, j+1:nbPoints) = 1/(2*pi)* ...
                    XdotY(j+1:nbPoints)./norm_XY_square(j+1:nbPoints).*sigma(j+1:nbPoints);
                
                % Diagonal terms
                NPOpeKerMat(j,j) = -1/(4*pi)* ...
                    avec(:,j)'*normal(:,j)/(tvec_norm(j)^2)*sigma(j);
            end
        end    
        
        
    end 
end