classdef beamformed_data < handle
    %BEAMFORMED_DATA   beamformed_data definition. Children of HANDLE class
    %
    %   See also PULSE, BEAM, PHANTOM, PROBE
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %   $Date: 2017/03/08 $
    
    %% compulsory properties
    properties  (SetAccess = public)
        scan                       % SCAN class
        data                       % data
    end
    
    %% optional properties
    properties  (SetAccess = public)
        phantom                    % PHANTOM class [optional]
        wave                       % WAVE class [optional]
        probe                      % PROBE class [optional]
        pulse                      % PULSE class [optional]
    end
    
    %% dependent properties
    properties  (Dependent)
    end
    
    
    %% constructor
    methods (Access = public)
        function h=beamformed_data()
            %BEAMFORMED_DATA   Constructor of beamformed_data class
            %
            %   Syntax:
            %   h = beamformed_data()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE
        end
    end
    
    %% copy
    methods (Access = public)
        function copy(h,object)
            %COPY    Copy the values from another BEAMFORMED_DATA class
            %
            %   Syntax:
            %   COPY(object)
            %       object       Instance of a BEAMFORMED_DATA class
            %
            %   See also SCAN, WAVE, SOURCE
            assert(isa(object,class(h)),'Class of the input object is not identical');
            
            % we copy all non-dependent public properties
            list_properties=properties(object);
            for n=1:numel(list_properties)
                property_name=list_properties{n};
                mp = findprop(h,property_name);
                if strcmp(mp.GetAccess,'public')&&~mp.Dependent
                    eval(sprintf('h.%s = object.%s',property_name,property_name));
                end
            end
        end
    end
    
    %% plot methods
    methods
        function figure_handle=plot(h,figure_handle_in,in_title,dynamic_range)
           
            if nargin>1 && ~isempty(figure_handle_in)
                figure_handle=figure(figure_handle_in);
            else
                figure_handle=figure();
            end
            if nargin<3 
                in_title='Beamformed data';
            end
            if nargin<4 
                dynamic_range=60;
            end
            
            % convert to intensity values
            envelope=abs(h.data);
            envelope_dB=20*log10(envelope./max(envelope(:)));
            
            switch class(h.scan)
                case 'huff.linear_scan'
                    if size(envelope_dB,3)>1
                        envelope_dB=max(envelope_dB,-dynamic_range);
                        envelope_dB=(envelope_dB+dynamic_range)/dynamic_range;
                        implay(reshape(envelope_dB,[h.scan.N_z_axis h.scan.N_x_axis 1 size(envelope_dB,3)]),10);                        
                    else
                        x_matrix=reshape(h.scan.x,[h.scan.N_z_axis h.scan.N_x_axis]);
                        z_matrix=reshape(h.scan.z,[h.scan.N_z_axis h.scan.N_x_axis ]);
                        pcolor(x_matrix*1e3,z_matrix*1e3,reshape(envelope_dB,[h.scan.N_z_axis h.scan.N_x_axis]));
                        %imagesc(h.scan.x_axis*1e3,h.scan.z_axis*1e3,reshape(envelope_dB,[h.scan.N_z_axis h.scan.N_x_axis]));
                        shading flat;
                        set(gca,'fontsize',14);
                        set(gca,'YDir','reverse');
                        axis tight equal; 
                        colorbar; 
                        colormap gray;
                        xlabel('x[mm]'); ylabel('z[mm]');
                        caxis([-dynamic_range 0]);
                        title(in_title);
                    end
                case 'huff.sector_scan'
                    if size(envelope_dB,3)>1
                        envelope_dB=max(envelope_dB,-dynamic_range);
                        envelope_dB=(envelope_dB+dynamic_range)/dynamic_range;
                        implay(reshape(envelope_dB,[h.scan.N_azimuth_axis h.scan.N_depth_axis 1 size(envelope_dB,3)]),10);                        
                    else
                        x_matrix=reshape(h.scan.x,[h.scan.N_depth_axis h.scan.N_azimuth_axis]);
                        z_matrix=reshape(h.scan.z,[h.scan.N_depth_axis h.scan.N_azimuth_axis ]);
                        pcolor(x_matrix*1e3,z_matrix*1e3,reshape(envelope_dB,[h.scan.N_depth_axis h.scan.N_azimuth_axis]));
                        shading flat;
                        set(gca,'fontsize',14); 
                        set(gca,'YDir','reverse');
                        axis tight equal; 
                        colorbar; 
                        colormap gray;
                        xlabel('x[mm]'); ylabel('z[mm]');
                        caxis([-dynamic_range 0]);
                        title(in_title);
                    end
                otherwise
                    error(sprintf('Dont know how to plot on a %s yet. Sorry!',class(h.scan)));
            end
        end
    end
    
    %% set methods
    methods
        function h=set.phantom(h,in_phantom)
            if ~isempty(in_phantom)
                assert(isa(in_phantom,'huff.phantom'), 'The input is not a PHANTOM class. Check HELP PHANTOM.');
                h.phantom=in_phantom;
            end
        end
        function h=set.pulse(h,in_pulse)
            if ~isempty(in_pulse)
                assert(isa(in_pulse,'huff.pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
                h.pulse=in_pulse;
            end
        end
        function h=set.probe(h,in_probe)
            if ~isempty(in_probe)
                assert(isa(in_probe,'huff.probe'), 'The input is not a PROBE class. Check HELP PROBE.');
                h.probe=in_probe;
            end
        end
        function h=set.wave(h,in_wave)
            if ~isempty(in_wave)
                assert(isa(in_wave,'huff.wave'), 'The input is not a WAVE class. Check HELP WAVE.');
                h.wave=in_wave;
            end
        end
        function h=set.scan(h,in_scan)
            if ~isempty(in_scan)            
                assert(isa(in_scan,'huff.scan'), 'The input is not a SCAN class. Check HELP SCAN.');
                h.scan=in_scan;
            end
        end
        function h=set.data(h,in_data)
            % some checking is due here
            
            h.data=in_data;
        end
    end
    
    %% get methods
    methods
        
    end
end