function vectorEst = vectorDoppler(R1,p)
%Gives you the x and z velocity components, p.angles are the steering
%angles with respect to the transducer normal. Assuming receive direction
%normal to transducer.

sizeR1 = size( squeeze(R1) );

if sizeR1(2) ~= 1
    vectorEst = single(zeros(size(R1,1), size(R1,2),2));
    temp = vectorEst;
    temp = reshape(temp, [size(temp,1)*size(temp,2),size(temp,3)]);
    fMat = (reshape(angle(R1),[size(R1,1)*size(R1,2),size(R1,3)])*p.PRF/(2*pi) ).';

    %aMat = [-sin(p.angles/2); cos(p.angles/2)]; % Works as an approximation to the truth for small steering angles
    aMat = [-sin(p.angles) (1+cos(p.angles))]./2; 

    pseudoInv = pinv(aMat).';
    for kk = 1:size(fMat,2)

       temp(kk,:) = pseudoInv*fMat(:,kk);
    end


    vectorEst = reshape(temp,[size(vectorEst,1), size(vectorEst,2), size(vectorEst,3)]);

else
    fMat = ( angle(R1)*p.PRF/(2*pi) ).';

    %aMat = [-sin((p.angles/2).'); cos((p.angles/2).')]; % Works for small steering angles
     aMat = [-sin(p.angles).'; (1+cos(p.angles)).']./2;

    pseudoInv = pinv(aMat).';
    vectorEst = pseudoInv*fMat';
 
    
end

end