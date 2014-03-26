function v = wfPhononCoordinates(wf, uIr, uRam, nIr, nRam)

% Projects the wavefunction wf into phonon coordinates (uIr, uRam)
% (Write projection equation and references here)

pkg load miscellaneous

v = 0.0;
for e1 = 1:3
  for e2 = 1:3
    partialSum = 0.0;
    for ir = 0:nIr
      for ram = 0:nRam
	partialSum += wf(stateLabel(e1, e2, ir, ram, nIr)) \
	              .* (1.0/sqrt(2^(ir+ram) * factorial(ir) * \
				  factorial(ram) * pi)) \
                      .* exp(- (uIr.^2 + uRam.^2) ./ 2.0) \
	              .* hermitepoly(ir, uIr) .* hermitepoly(ram, uRam);
      endfor
    endfor
    v += (partialSum.^2);
  endfor
endfor

end
