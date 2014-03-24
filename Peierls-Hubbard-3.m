%% Computing the eigenvalues and some more features 
%  of a 3-sites Peierls-Hubbard quantum model

% (Write hamiltonian and references here)

% Relevant parameters
band_energies = [0.5 -0.5 0.5];
nn_hopping = 0.5;
on_site_repulsion = 7.0;
ir_energy = 612.4 / 8056;
ram_energy = 500 / 8056;
ram_shift = 4.0 / 3.0;
e_ir_coupling = 0.13;
e_ram_coupling = 0.0;
n_ir = 20;
n_ram = 20;

size = 9 * (n_ir + 1) * (n_ram + 1);

% To build the hamiltonian we store the non-zero values as triplets (row, col, val)
rows = [];
cols = [];
vals = [];

for e1 = 1:3
  for e2 = 1:3
    for ir = 0:n_ir
      for ram = 0:n_ram
	
	% Band energies
	rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	vals(end + 1) = band_energies(e1) + band_energies(e2);

	% On-site Coulomb repulsion
	if (e1 == e2)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  vals(end + 1) = on_site_repulsion;
	endif

	% Nearest-neighbor hopping
	if (e1 != 3)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1 + 1, e2, ir, ram, n_ir);
	  vals(end + 1) = nn_hopping;
	endif
	if (e1 != 1)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1 - 1, e2, ir, ram, n_ir);
	  vals(end + 1) = nn_hopping;
	endif
	if (e2 != 3)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2 + 1, ir, ram, n_ir);
	  vals(end + 1) = nn_hopping;
	endif
	if (e2 != 1)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2 - 1, ir, ram, n_ir);
	  vals(end + 1) = nn_hopping;
	endif

	% Phonons' energies
	rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	vals(end + 1) = ir * ir_energy;

	rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	cols(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	vals(end + 1) = ram * ram_energy;

	% Electron - infrared phonons coupling
	n = e1 + e2 - 4;
	if (ir != n_ir)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2, ir + 1, ram, n_ir);
	  vals(end + 1) = n * e_ir_coupling * sqrt(ir + 1);
	endif
	if (ir != 0)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2, ir - 1, ram, n_ir);
	  vals(end + 1) = n * e_ir_coupling * sqrt(ir);
	endif
	
	% Electron - Raman phonons coupling
	n = abs(e1 - 2) + abs(e2 - 2) - ram_shift;
	if (ram != n_ram)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram + 1, n_ir);
	  vals(end + 1) = n * e_ram_coupling * sqrt(ram + 1);
	endif
	if (ram != 0)
	  rows(end + 1) = stateLabel(e1, e2, ir, ram, n_ir);
	  cols(end + 1) = stateLabel(e1, e2, ir, ram - 1, n_ir);
	  vals(end + 1) = n * e_ram_coupling * sqrt(ram);
	endif
      endfor
    endfor
  endfor
endfor

h = spconvert([rows; cols; vals]');

