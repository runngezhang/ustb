classdef scpw < us_dataset
%SCPW    Symmetric coherent plane wave (sta) dataset.
%
%   See also US_DATASET, SCPW.SCPW

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/03/05 $
    
    properties (SetAccess = public)
        tx_angle       % vector containing the transmit angles (rad)
        rx_angle       % vector containing the receive angles (rad)
    end
    
    %% constructor
    methods (Access = public)       
        function h = scpw(name,input_format,input_c0,input_tx_angle,input_rx_angle,input_time,input_data,input_geom,input_modulation_frequency)
            %CPW    Constructor of the sta class.
            %
            %   Syntax:
            %   CPW(name,format,c0,angle,time,data,geom,modulation_frequency) 
            %       name                    Name of the dataset
            %       format                  Format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %       c0                      Reference speed of sound (m/s)
            %       tx_angle                Transmit Plane wave angle vector (rad)
            %       rx_angle                Receive Plane wave angle vector (rad)
            %       time                    Time vector (s)
            %       data                    Numerical data [time_samples, channels, firings, frames]
            %       geom                    Probe geometry [x, y, z] (m)
            %       modulation_frequency    Modulation frequency (Hz) - required only for IQ format
            %
            %   See also CPW, US_DATASET
        
            % we allow object instantiation without parameters
            h@us_dataset(); 
            
            if exist('name') 
                % required data
                h.name = name;
                h.format = input_format;
                h.c0 = input_c0;
                h.time = input_time;
                h.data = input_data;
                h.tx_angle = input_tx_angle;
                h.rx_angle = input_rx_angle;
                h.geom = input_geom;

                if(h.format==E.signal_format.IQ) 
                    h.modulation_frequency=input_modulation_frequency;
                end

                % checks
                assert(size(h.data,2)==length(h.rx_angle),'The number of receive angles do not match data format');
                assert(size(h.data,3)==length(h.tx_angle),'The number of transmit angles do not match data format');
            end
        end
    end
        
    %% image reconstruction
    methods  (Access = public)       
         function elapsed_time = image_reconstruction(h,recons,implem)
            %IMAGE_RECONSTRUCTION    Method for image reconstruction of Coherent plane wave datasets
            %
            %   Syntax:
            %   IMAGE_RECONSTRUCTION(transmit_beam,receive_beam,scan,implementation)
            %       reconstruction          Class containing the reconstruction specification
            %       implementation          Enumeration specifying the algorithm 
            %
            %   See also RECONSTRUCTION, CPW
             
            % default implementation
            if ~exist('implem') implem=E.implementation.mex; end
            
            % initial time -> performance test
            tic; 
            
            % loop over orientations
            total_data=zeros(size(recons.scan.x_matrix,1),size(recons.scan.x_matrix,2),length(recons.orientation),size(h.data,4));
            for o=1:length(recons.orientation)
                ef=0.3;                                                          % half fractional bandwidth (hardcoded for the moment)
                % transmit apodization
                xT=recons.scan.x*ones(1,length(h.tx_angle))-recons.scan.z*tan(h.tx_angle.'); % position of equivalent receive element -> Alfonso's equation 
                pf=-erf((abs(sin(h.tx_angle))*2*recons.orientation.transmit_beam.f_number-1.0)/ef).';
                Pf=ones(recons.scan.pixels,1)*pf;                                % Montaldo's wideband extension
                h.transmit_apodization = Pf.* recons.calculate_apodization(recons.orientation(o).transmit_beam,xT,0.*xT);
                % receive apodization
                pf=-erf((abs(sin(h.rx_angle))*2*recons.orientation.transmit_beam.f_number-1.0)/ef).';
                Pf=ones(recons.scan.pixels,1)*pf;                                % Montaldo's wideband extension
                xR=recons.scan.x*ones(1,length(h.rx_angle))-recons.scan.z*tan(h.rx_angle.'); % position of equivalent receive element -> Alfonso's equation 
                h.receive_apodization = recons.calculate_apodization(recons.orientation(o).receive_beam,xR,0.*xR);
                % launch selected implementation
                temporal_data=h.launch_implementation(recons,implem);
                
                % reshape matrix to include orientation dimensions
                total_data(:,:,o,:)=reshape(temporal_data,[size(recons.scan.x_matrix,1) size(recons.scan.x_matrix,2) 1 size(temporal_data,2)]);
            end
                        
            % copy data to reconstruction
            recons.data=total_data;
            recons.format=h.format;
            if(h.format==E.signal_format.RF)
                [recons.central_frequency, recons.bandwidth]=tools.estimate_frequency(h.time,h.data);
            end
                        
            % elapsed time for performance test
            elapsed_time=toc; 
        end
    end

    %% set methods
    methods  
        function set.tx_angle(h,input_angle)
            assert(size(input_angle,1)>=size(input_angle,2), 'The angle must be a column vector');
            h.tx_angle=input_angle;
        end
        function set.rx_angle(h,input_angle)
            assert(size(input_angle,1)>=size(input_angle,2), 'The angle must be a column vector');
            h.rx_angle=input_angle;
        end
    end
    
end
