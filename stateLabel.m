function label = stateLabel(e1, e2, ir, ram, nIr)

label = e1 + (3 * (e2 - 1)) + (9 * ir) + (9 * ram * (nIr +1));

end
