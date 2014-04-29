function [bf X Z] = planewave_beamforming2(ch_data,p,x,z,apod)
% bf = planewave_beamforming2(ch_data,p,x,z,apod)
% 
% Beamforming method for plane wave emissions. The parameter structure 'p'
% should have the follwoing fields:
% c         - Speed f sound [m/s]
% dx        - Array pitch
% t0        - Time of first sample
% tx_angle  - Transmit angle in radians
% rx_angle  - Receive angle in radians
% FN        - F number for receive aperture.
% fs_in     - Input Sampling frequency
% fs_out    - Output Sampling frequency
% f_demod   - IQ demodulation frequency. Needed if data is complex
% 
% Input:
% rf_ch     - Channel data, RF or IQ
% p         - Parameter structure
% x         - x coordinates for beamforming area. If receive steering is
% applied this is the x coordinate of the starting point if each receive
% line. The offset is computed
% z         - z coordinates for beamforming area. If receive steering is
% applied this is the radial coordinates.
% apod      - [Optional] 0 = rectangular apodization, 1 = hamming
% apodization
%
% Output
% bf        - Beamformed data.
% X         - x-grid coordinates for the beamformed data
% Z         - z-grid coordinates for the beamformed data
% 