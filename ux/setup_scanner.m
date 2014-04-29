function success = setup_scanner(ult,sc_setup)

success = 1;
for kk=1:length(sc_setup.param_ids)
    if ~ult.set_int_param(sc_setup.param_ids(kk),sc_setup.param_values(kk))
        fprintf('Could not set parameter with ID : %d and value : %d \n',sc_setup.param_ids(kk),sc_setup.param_values(kk));
        success = 0;
    end
end
