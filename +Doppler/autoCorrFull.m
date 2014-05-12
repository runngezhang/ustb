function R1 = autoCorrFull(data_path,fName,p,nx,nz)
%%% Currently hacked to fit with our data...should be fixed

import Doppler.*

fullFile = matfile(fullfile(data_path,fName));

R1 = single(zeros(p.R,p.B,length(p.angles)));
for kk = 1:length(p.angles)
    iqTemp = squeeze(fullFile.IQ_full(:,:,kk:15:end));
    iqTemp = permute(iqTemp(:,:,p.Nstart:p.Nstart+p.packetSize-1),[3,1,2]);
    [iqhp vThresh] = hp(iqTemp,p);
    clear iqTemp
    R1(:,:,kk) = squeeze(mean(conj(iqhp(1:end-1,:,:)).*iqhp(2:end,:,:))); % R1 estimate
    clear iqhp
    
    if nargin>3
        myMask = (1/(nz*nx))*ones(nz,nx);
        R1(:,:,kk) = filter2(myMask,R1(:,:,kk));
    end
    
end

end
