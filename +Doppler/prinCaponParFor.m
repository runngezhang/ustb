function [P,lambvec,lambvec_alt] = prinCaponParFor(Nf,ROI,s)
nSig = s.nSig;
diagonalLoading = s.diagonalLoading;
setSig = s.setSig;
fbavg = s.fbavg;
clutterPresence = s.clutter;
estSigMethod = s.estSigMethod;

CC = ROI(:,:,:,:);
[K, R, B, F] = size(CC);
P = zeros(Nf,F);

lambvec = zeros(1,F);
matlabpool(10)
parfor l = 1:F %for all frames
    a = 0;
    frame = CC(:,:,:,l);    
    Rxx = (zeros(K,K,R));

    for j = 1:B %for all beams
        plane = squeeze(frame(:,:,j));
        for k = 1:R %for all ranges
            x = plane(:,k);
            Rxx(:,:,k+a) = x*x';
        end
        a = a+R;
    end

    Rxxm = mean(Rxx,3);
    
    if fbavg
        J  = eye(K); J = fliplr(J); 
        Rxxm  = 1/2*(Rxxm + J*Rxxm.'*J);  %As done in adaptive beamforming  
    else
        Rxxm = 0.5*(Rxxm' + Rxxm); % Force real diagonal, and Hermetian symmetry  (Equal to backward forward averaging?)
    end
    
   
    if diagonalLoading %Diagonal loading: add small positive number to diagonal. Proportional to power in frame
        RxxDiag = Rxxm + (1/(1000*size(Rxxm,1)))*sum(diag(Rxxm))*eye(size(Rxxm)); 
    else
        RxxDiag = Rxxm;
    end
    
    
    [Evecs,Evals] = eig(RxxDiag);
    [Evals,Index] = sort(diag(Evals),'descend');
    Evecs         = Evecs(:,Index);        
    R_            = Evecs*diag(1./Evals)*Evecs'; %Estimate of inverse covariance matrix, using all components, noise and signal subspace
    
    f = (freqspace(Nf,'whole')-1)*0.5; 
    k = (0:K-1);
    Pff = zeros(length(f),1);
    q = complex(0,1);
    %Need to estimate how many signal components are in each frame
    
    maxsum = 0;
    Lambda_no = 0;
    sEProj = 0;
   
   
    %% Method one for finding number of signal components in each frame
     
    for i = 1:length(f);
         e = 2*f(i)*pi*k;
         Ef = exp(q*e).'; 
         maxsum  = maxsum + sum(abs(Evecs(:,1:K)'*Ef).^2./(Evals(1:K)));
    end
    
    if diagonalLoading
        threshP = 0.005;
    else
        threshP = 0.0001;
    end
    
    while sEProj < (threshP*maxsum) %Find suitable threshold value to insert here
        Lambda_no   = Lambda_no+1;
        sEProj = 0;
        for i = 1:length(f);
             e = 2*f(i)*pi*k;
             Ef = exp(q*e).' ;  
             sEProj  = sEProj + sum(abs(Evecs(:,1:Lambda_no)'*Ef).^2./(Evals(1:Lambda_no)));
        end
         %clear e Ef;
    end
 

    %% Method two for finding number of signal components in each frame
    maxSumEvals = sum(Evals(1:K));
    if clutterPresence
        Threshold = 0.995; %Some percentage of the signal energy is enough
    else
        Threshold = 0.5;
    end
    sumEvals = 0;
    for ix = 1:K
        sumEvals = sumEvals + Evals(ix);
        if sumEvals>=Threshold*maxSumEvals
            break
        end
    end
    Lambda_no_alt = ix;
    %% Force number of signal components?
    if setSig == 1
        Lambda_no = nSig;  
        Lambda_no_alt = nSig;
    else
         if clutterPresence ==1;
            if Lambda_no<2 || Lambda_no_alt<2 
                Lambda_no = 2;
                Lambda_no_alt = 2;
                disp('Forcing at least two signal components-only suitable for non filtered data!')
            end
         else
             if Lambda_no>4 || Lambda_no_alt>4 
                 Lambda_no = 4;
                 Lambda_no_alt = 4;
             end
         end
    end
    
    %% Calculate power estimate through estimate of the signal amplitude
    lambvec(l) = Lambda_no; %To look at how the number of estimated signal components vary with every frame
    lambvec_alt(l) = Lambda_no_alt;
    
    switch estSigMethod
        case 'energyPercentage'
            EvecsSig=Evecs(:,1:Lambda_no_alt); %Use either set or estimated nr of signal components in projection of weight vector
            disp(sprintf('Using %d percent of the signal energy; method two for finding number of signal components',Threshold*100))
        case 'energyProjection'
            disp(sprintf('Using projection technique with treshold %d; method one for finding number of signal components',threshP))
            EvecsSig=Evecs(:,1:Lambda_no);
    end
    
    for i = 1:length(f);
        e = 2*f(i)*pi*k;
        Ef = exp(q*e).'; 
        h = (R_*Ef)/( Ef'*R_*Ef) ;  
        hSig = EvecsSig*EvecsSig'*h; %Project optimum weights vector into signal subspace
        Pff(i) = hSig'*Rxxm*hSig; %Use optimal weights on original covariance matrix estimate without loading og subspace proj.
    end

    P(:,l) = Pff; %Power spectrum for all frames
end

matlabpool close
% figure(), plot(lambvec,'bo');
% hold on; plot(lambvec_alt,'ro');
% hold off;
end
