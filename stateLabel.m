function label = stateLabel(e1, e2, ir, ram, n_ir)

label = e1 + (3 * (e2 - 1)) + (9 * ir) + (9 * ram * (n_ir +1));

end
