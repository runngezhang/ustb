function pyusrecons = create_pyusrecons(recons)

tx_beam = recons.orientation.transmit_beam;
rx_beam = recons.orientation.receive_beam;

pyusrecons = py.pyus.Interop.USTBInterop.Reconstruction();
pyusrecons.name = recons.name;
pyusrecons.scan = py.pyus.Interop.USTBInterop.LinearScan();
pyusrecons.scan.x_axis = recons.scan.x_axis';
pyusrecons.scan.z_axis = recons.scan.z_axis';
pyusrecons.scan.update_attributes();
pyusrecons.orientation{1}.transmit_beam = py.pyus.Interop.USTBInterop.Beam(tx_beam.f_number, int32(tx_beam.apodization));  
pyusrecons.orientation{1}.receive_beam = py.pyus.Interop.USTBInterop.Beam(rx_beam.f_number, int32(rx_beam.apodization));

end

