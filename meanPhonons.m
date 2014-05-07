function [meanIr, meanRam, stddevIr, stddevRam] = meanPhonons(vecs, nIr, nRam)

%% Calculates the mean infrared and raman phonons as well as
%  their standard deviations for the linearPeierlsHubbard model

assert (size(vecs)(1) == 9 * (nIr + 1) * (nRam + 1), "Eigenvector matrix's size doesn't correspond with the provided number of phonons.")

% Number of eigenvectors
n = size(vecs)(2);

meanIr = zeros(1,n);
meanRam = zeros(1,n);
stddevIr = zeros(1,n);
stddevRam = zeros(1,n);

for i = 1:n
  sqrIr = 0;
  sqrRam = 0;
  for e1 = 1:3
    for e2 = 1:3
      for ir = 0:nIr
	for ram = 0:nRam
	  meanIr(i) += ir * vecs(stateLabel(e1,e2,ir,ram,nIr), i)^2;
	  meanRam(i) += ram * vecs(stateLabel(e1,e2,ir,ram,nIr), i)^2;
	  sqrIr += (ir^2) * vecs(stateLabel(e1,e2,ir,ram,nIr),i)^2;
	  sqrRam += (ram^2) * vecs(stateLabel(e1,e2,ir,ram,nIr),i)^2;
	endfor
      endfor
    endfor
  endfor
  stddevIr(i) = sqrt(sqrIr - meanIr(i)^2);
  stddevRam(i) = sqrt(sqrRam - meanRam(i)^2);
endfor
