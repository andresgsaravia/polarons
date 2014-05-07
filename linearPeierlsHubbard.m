function [h, eigvals, eigvecs] = linearPeierlsHubbard(bandEnergies, 
						      nnHopping, onSiteRepulsion,
						      irEnergy, ramEnergy, ramShift, eIrCoupling,
						      eRamCoupling,
						      nIr, nRam, k)

%% Computing the eigenvalues and some more features 
%  of a 3-sites Peierls-Hubbard quantum model

% (Write hamiltonian and references here)

% To build the hamiltonian we store the non-zero values as triplets (row, col, val)
rows = [];
cols = [];
vals = [];

for e1 = 1:3
  for e2 = 1:3
    for ir = 0:nIr
      for ram = 0:nRam
	
	% Band energies
	rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	vals(end + 1) = bandEnergies(e1) + bandEnergies(e2);

	% On-site Coulomb repulsion
	if (e1 == e2)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  vals(end + 1) = onSiteRepulsion;
	endif

	% Nearest-neighbor hopping
	if (e1 != 3)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1 + 1, e2, ir, ram, nIr);
	  vals(end + 1) = nnHopping;
	endif
	if (e1 != 1)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1 - 1, e2, ir, ram, nIr);
	  vals(end + 1) = nnHopping;
	endif
	if (e2 != 3)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2 + 1, ir, ram, nIr);
	  vals(end + 1) = nnHopping;
	endif
	if (e2 != 1)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2 - 1, ir, ram, nIr);
	  vals(end + 1) = nnHopping;
	endif

	% Phonons' energies
	rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	vals(end + 1) = ir * irEnergy;

	rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	vals(end + 1) = ram * ramEnergy;

	% Electron - infrared phonons coupling
	n = e1 + e2 - 4;
	if (ir != nIr)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2, ir + 1, ram, nIr);
	  vals(end + 1) = n * eIrCoupling * sqrt(ir + 1);
	endif
	if (ir != 0)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2, ir - 1, ram, nIr);
	  vals(end + 1) = n * eIrCoupling * sqrt(ir);
	endif
	
	% Electron - Raman phonons coupling
	n = abs(e1 - 2) + abs(e2 - 2) - ramShift;
	if (ram != nRam)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram + 1, nIr);
	  vals(end + 1) = n * eRamCoupling * sqrt(ram + 1);
	endif
	if (ram != 0)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, nIr);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram - 1, nIr);
	  vals(end + 1) = n * eRamCoupling * sqrt(ram);
	endif
      endfor
    endfor
  endfor
endfor

% Building the hamiltonian
h = spconvert([rows; cols; vals]');
[eigvecs, eigvals] = eigs(h, k, "sa");

end
