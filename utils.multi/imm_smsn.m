function out = imm_smsn(y, model)

addpath('../dens');
addpath('../utiles');

adjoint = @(A) det(A)*inv(A);
deriv_der = @(A,B,C) det(A)*sum(B * transpose(C), 'all');

family = model.class;

if((family ~= "t") && (family ~= "Skew.t") && (family ~= "Skew.cn") && (family ~= "Skew.slash") && (family ~= "Skew.normal") && (family ~= "Normal"))
    error("Family "+ family +" not recognized.")
end

if (size(y,2) < 2)
    error("The dimension of y (p) is: " + num2str(size(y,2)) +". We need p >= 2.")
end

if (model.uni_Gama)
    error("Sorry. The information matrix cannot be calculated when the uni.Gama was used!")
end

[n,p] = size(y);
g = length(model.pii);

% all model parameters are cell

for i = 1 : length(model.Sigma)
    model.Sigma{i} = model.Sigma{i} * model.Sigma{i};
end

mu = model.mu;
Sigma = model.Sigma;
lambda = model.shape;
pii = model.pii;
nu = model.nu;

if (class(model) == "t")
    soma = 0;
    soma2 = 0;
    
    I_Phi = @(w,Ai,di,nu) ((2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*tcdf(((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu);
    I_phi = @(w,Ai,di,nu) ((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2);
    
    for i = 1 : n
        S = [];
        dPsi_dnu = 0;
        yi = y(i,:);
        for j = 1 : g
            Dr = matrix_sqrt(Sigma{j});
            Dr_inv = inv(Dr);
            d_sig = det(Sigma{j});
            
            Ai = lambda{j}*Dr.inv*transpose(y(i,:) - mu{j});
            di = mahalanobis(yi, mu{j}, Sigma{j});
            
            
            %derivative
            dir_dmu = -2*(Dr_inv*Dr_inv)*trnspose(y(i,:) - mu{j});
            dAir_dmu = -Dr_inv*transpose(lambda{j});
            
            dPsi_dmu = ((2*d_sig^(-1/2))/(2*pi)^(p/2))*(dAir_dmu * I_phi((p+1)/2, Ai, di, nu) - (1/2)*dir_dmu*I_Phi((p/2)+1, Ai, di, nu));
            
            % for sigma elements
            l = 1; m = 1;
            for k = 1 : ((p+1)*p/2)
                Vis = false;
                D = zeros(p,p);
                D(l,m) = 1; D(m,l) = 1;
                
                ddet_ds = -(1/det(Dr)^2)*deriv_der(Dr,Dr_inv,D);
                dir_ds = - (y(i,:) - mu{j})*Dr.inv*(D*Dr.inv + Dr.inv*D)*Dr.inv*transpose(y(i,:) - mu{j});
                dAir_ds = - lambda{j}*Dr_inv*D*Dr.inv*transpose(y(i,:) - mu{j});
                
                dPsi_dsigma = (2/(2*pi)^(p/2))*(ddet_ds*I_Phi(p/2, Ai, di, nu) - (1/2)*dir_ds*d_sig^(-1/2)*I_Phi(p/2+1, Ai, di, nu) ...
                    + d_sig^(-1/2)*dAir_ds*I_phi((p+1)/2, Ai, di, nu));
                Ssigma(k) = (pii(j) / d_mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi_dsigma;
                
                if ((( l*m - p*floor((l*m)/p)) == 0) && (l ~= m))
                    l = l + 1;
                    m = l;
                    Vis = TRUE;
                end
                
                if(~Vis)
                    m = m+1;
                end
            end
            
            
            ui = random('gamma', nu/2, 2/nu,1,1e4);
            resto = mean(ui.^(p/2).*log(ui).*exp(-ui.*di/2).*normcdf(ui.^(1/2)*Ai));
            dPsi_dnu = dPsi_dnu + pii(j)*((d_sig^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-psi(nu/2))*I_Phi(p/2, Ai, di, nu) - I_Phi((p+2)/2, Ai, di, nu) + resto);
            
            Simu = (pii(j)/ d_mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi_dmu;
            
            
            
            S = [S, Simu, Ssigma];
        end
        Sinu = (1/d_mixedmvST(yi, pii, mu, Sigma, lambda, nu))*dPsi_dnu;
        if (g > 1)
            for j = 1 : (g-1)
                Sipi(j) = (1/d_mixedmvST(yi, pii, mu, Sigma, lambda, nu))*(dmvt_ls(yi, mu{j}, Sigma{j}, lambda{j}, nu) - dmvt_ls(yi, mu{g}, Sigma{g}, lambda{g}, nu));
                S = [S, Sipi, Sinu];
            end
        end
        if(g == 1)
            S = [S, Sinu];
        end
        soma = soma + transpose(S)*S;
        
    end
    
        NAME = []; piiN = [];
    for(i = 1 : g){ %I
        SigmaN = []; muN = []; shapeN = [];
        for (k = 1:p){
           muN = c(muN,paste("mu",i,"_",k,sep=""))
        end
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    } % I
    if( g==1) NAME <- c(NAME,"nu")
    else{    
       for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN,"nu")
    }

    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME
    
end

rmpath('../dens');
rmpath('../utiles');

