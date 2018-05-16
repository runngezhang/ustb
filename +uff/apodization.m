classdef apodization < uff
    %APODIZATION   UFF class to hold apodization data
    %   APODIZATION contains data to define transmit, receive & synthetic 
    %   beams. Different parameters are needed depending on the use.
    %
    %   Properties:
    %         probe               % UFF.PROBE class (needed for transmit & receive apodization)
    %         focus               % UFF.SCAN class (needed for transmit, receive & synthetic apodization)       
    %         sequence            % collection of UFF.WAVE classes (needed for synthetic apodizaton)
    %         
    %         window              % UFF.WINDOW class, default uff.window.noen
    %         f_number            % F-number [Fx Fy] [unitless unitless] 
    %         M                   % Number of elements [Mx My] in case f_number=0
    %         
    %         origo               % POINT class
    %         tilt                % tilt angle [azimuth elevation] [rad rad] 
    %
    %   Example:
    %         apo = uff.apodization();
    %
    %   See also UFF.CHANNEL_DATA, UFF.BEAMFORMED_DATA, UFF.SCAN
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %   $Last updated: 2017/06/09$

    %% public properties
    properties  (Access = public)
        probe                           % UFF.PROBE class (needed for transmit & receive apodization)
        focus                           % UFF.SCAN class (needed for transmit, receive & synthetic apodization)       
        sequence                        % collection of UFF.WAVE classes (needed for synthetic apodizaton)

        f_number  = [1 1]               % F-number [Fx Fy] [unitless unitless] 
        window    = uff.window.none     % UFF.WINDOW class, default uff.window.none
        MLA       = 1                   % number of multi-line acquisitions, only valid for uff.window.scanline
        MLA_overlap = 0                 % number of multi-line acquisitions, only valid for uff.window.scanline
        
        origo     = uff.point('xyz',[0, 0, -Inf]);  % POINT class
        tilt      = [0 0]                           % tilt angle [azimuth elevation] [rad rad] 
        minimum_aperture = [1e-3 1e-3]  % minimum aperture size in the [x y] direction
        maximum_aperture = [10 10]      % maximum aperture size in the [x y] direction
    end
    
    %% optional properties
    properties  (Access = public)
        apodization_vector              % apodization vector to override the dynamic calculation of apodization 
    end
    
    %% dependent properties
    properties  (Dependent)   
        data                        % apodization data 
        propagation_distance        % distance from the origin to the pixel
        N_elements                  % number of elements (real or synthetic)
    end
    
    %% private & Hidden
    properties  (Dependent, Hidden, Access = private)   
            origin                      % beam origin
    end
    
    %% private properties
    properties (Access = private)
        data_backup
    end
    
    %% constructor
    methods (Access = public)
        function h=apodization(varargin)
            h = h@uff(varargin{:});
        end
    end
    
    %% set methods
    methods  
        function h=set.origo(h,in_source)
            assert(isa(in_source,'uff.point'), 'The input origo is not a POINT class. Check HELP POINT');
            h.origo=in_source;
        end
        function h=set.probe(h,in_probe)
            assert(isa(in_probe,'uff.probe'), 'The input probe is not a PROBE class. Check HELP PROBE');
            h.probe=in_probe;
        end
        function h=set.focus(h,in_scan)
            assert(isa(in_scan,'uff.scan'), 'The input focus is not a SCAN class. Check HELP SCAN');
            h.focus=in_scan;
        end
        function h=set.f_number(h,in_f_number)
            if(numel(in_f_number)==1) % we allow for escalar input
                 in_f_number=[in_f_number in_f_number];
             end 
            assert(size(in_f_number,1)==1&&size(in_f_number,2)==2,'The f-number must be a row vector [Fx Fy]');
            h.f_number=in_f_number;
        end        
        function h=set.tilt(h,in_tilt)
            if(numel(in_tilt)==1) % we allow for escalar input
                 in_tilt=[in_tilt 0];
             end 
            assert(size(in_tilt,1)==1,size(in_tilt,1)==2,'The tilt must be a row vector [azimuth elevation]');
            h.tilt=in_tilt;
        end  
        function h=set.window(h,in_window)
             assert(isa(in_window,'uff.window'),'The window input should be a WINDOW class. Check help WINDOW');
             h.window=in_window;
        end
        
        function h=set.minimum_aperture(h,in_ap)
            if(numel(in_ap)==1) % we allow for escalar input
                 in_ap=[in_ap in_ap];
             end 
            assert(size(in_ap,1)==1&&size(in_ap,2)==2,'The minimum aperture should be a row vector [Ax Ay]');
            h.minimum_aperture=in_ap;
        end     
        
        function h=set.maximum_aperture(h,in_ap)
            if(numel(in_ap)==1) % we allow for escalar input
                 in_ap=[in_ap in_ap];
             end 
            assert(size(in_ap,1)==1&&size(in_ap,2)==2,'The maximum aperture should be a row vector [Ax Ay]');
            h.maximum_aperture=in_ap;
        end     
    end
    
    methods
        % angular version
        function value=rectangular(h,ratio)
            value=double(ratio<=0.5);
        end
        function value=hanning(h,ratio)
            value=double(ratio<=0.5).*(0.5 + 0.5*cos(2*pi*ratio));
        end
        function value=hamming(h,ratio)
            value=double(ratio<=0.5).*(0.53836 + 0.46164*cos(2*pi*ratio));
        end
        function value=tukey(h,ratio, roll)            
            value=(ratio<=(1/2*(1-roll))) + (ratio>(1/2*(1-roll))).*(ratio<(1/2)).*0.5.*(1+cos(2*pi/roll*(ratio-roll/2-1/2)));
        end
    end
    
    %% get method
    methods
        function value=get.data(h)
            % check if we can skip calculation
            if h.check_hash()
                value = h.data_backup; 
                return;
            end
            
            % check if overridden
            if ~isempty(h.apodization_vector)
                value = h.apodization_vector;
                return;
            end
           
            % checking we have all we need
            if isempty(h.focus)
                h.focus=uff.scan('xyz',[0 0 0]);
            end
            if isempty(h.origo)
                h.origo=uff.point('xyz',[0 0 0]);
            end
            
            % NONE APODIZATION
            if(h.window==uff.window.none)
                h.data_backup=ones(h.focus.N_pixels,h.N_elements);
            elseif (h.window==uff.window.sta)
            % STA APODIZATION (just the element closest to origo)
                assert(numel(h.probe)>0,'The PROBE parameter must be set to use STA apodization.');
                dist=sqrt((h.probe.x-h.origo.x).^2+(h.probe.y-h.origo.y).^2+(h.probe.z-h.origo.z).^2);
                h.data_backup=ones(h.focus.N_pixels,1)*double(dist==min(dist(:)));
            elseif (h.window==uff.window.scanline)
            % SCALINE APODIZATION (MLA scanlines per wave)
                assert(numel(h.sequence)>0,'uff.apodization:Scanline','The SEQUENCE parameter must be set to use uff.window.scanline apodization.');
                if isa(h.focus,'uff.linear_scan')
                    N_waves=numel(h.sequence);
                    assert(N_waves==h.focus.N_x_axis/h.MLA,'The number of waves in the sequence does not match with the number of scanlines and set MLA.');
                    ACell=repmat({ones(h.MLA,1)},[1,h.focus.N_x_axis/h.MLA]);
                    if (h.MLA_overlap>0)
                        ABlock=filtfilt(ones(1,h.MLA_overlap+1)/(h.MLA_overlap+1),1,blkdiag(ACell{:}));
                    else
                        ABlock=blkdiag(ACell{:});
                    end
                    h.data_backup=kron(ABlock,ones(h.focus.N_z_axis,1));
                elseif isa(h.focus,'uff.sector_scan')  
                    N_waves=numel(h.sequence);
                    assert(N_waves==h.focus.N_azimuth_axis/h.MLA,'The number of waves in the sequence does not match with the number of scanlines and set MLA.');
                    ACell=repmat({ones(h.MLA,1)},[1,h.focus.N_azimuth_axis/h.MLA]);
                    if (h.MLA_overlap>0)
                        ABlock=filtfilt(ones(1,h.MLA_overlap+1)/(h.MLA_overlap+1),1,blkdiag(ACell{:}));
                    else
                        ABlock=blkdiag(ACell{:});
                    end
                    h.data_backup=kron(ABlock,ones(h.focus.N_depth_axis,1));
                else
                    error('uff.apodization:Scanline','The scan class does not support scanline based beamforming. This must be done manually, defining several scan and setting the apodization to uff.window.none.');
                end
            else

                % incidence ratio
                [ratio_theta ratio_phi]=h.incidence_ratios();
                
                % SWITCH
                switch(h.window)
                    % BOXCAR/FLAT/RECTANGULAR
                    case uff.window.boxcar
                        h.data_backup=h.rectangular(ratio_theta).*h.rectangular(ratio_phi);
                    % HANNING
                    case uff.window.hanning
                        h.data_backup=h.hanning(ratio_theta).*h.hanning(ratio_phi);
                    % HAMMING
                    case uff.window.hamming
                        h.data_backup=h.hamming(ratio_theta).*h.hamming(ratio_phi);
                    % TUKEY25        
                    case uff.window.tukey25
                        roll=0.25;
                        h.data_backup=h.tukey(ratio_theta,roll).*h.tukey(ratio_phi,roll);
                    % TUKEY50        
                    case uff.window.tukey50
                        roll=0.50;
                        h.data_backup=h.tukey(ratio_theta,roll).*h.tukey(ratio_phi,roll);
                    % TUKEY75        
                    case uff.window.tukey75
                        roll=0.75;
                        h.data_backup=h.tukey(ratio_theta,roll).*h.tukey(ratio_phi,roll);
                    % TUKEY80        
                    case uff.window.tukey80
                        roll=0.80;
                        h.data_backup=h.tukey_ang(ratio_theta,roll).*h.tukey_ang(ratio_phi,roll);
                    otherwise
                        error('Unknown apodization type!');
                end

    %             % edge smoothing
    %             if(size(apo,2)>(2*(beam.smoothing+1)))
    %                 % compute edge smoothing mask
    %                 mask=0.5-0.5*cos((0:beam.smoothing)/beam.smoothing*pi);  % vector mask
    %                 mask=(ones(h.focus.pixels,1)*mask(2:(beam.smoothing+1))); % matrix mask
    %                 
    %                 % modifying apodization edges
    %                 apo(:,1:beam.smoothing)=apo(:,1:beam.smoothing).*mask;
    %                 apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events)=apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events).*fliplr(mask);
    %             end
            end
            
            % normalize
            %h.data_backup=h.data_backup./sum(sum(h.data_backup,3),2);
            
            % pass a reference
            value=h.data_backup;
                
            % update hash
            h.save_hash();
        end
        
        function value=get.origin(h)
            % checking we have all we need
            %assert(numel(h.probe)>0,'The PROBE parameter is not set.');
            assert(h.focus.N_pixels>0,'The focus parameter is not set.');
            assert(numel(h.origo)>0,'The origo parameter is not set.');
            
            % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
            if isinf(h.origo.distance)
                origin(:,1)=h.focus.x+h.focus.z*tan(h.origo.azimuth)-h.focus.z*tan(h.tilt(1));
                origin(:,2)=h.focus.y+h.focus.z*tan(h.origo.elevation)-h.focus.z*tan(h.tilt(2));
            else
                origin(:,1)=h.origo.x-h.origo.z*(h.focus.x-h.origo.x)./(h.focus.z-h.origo.z)-h.focus.z*tan(h.tilt(1));
                origin(:,2)=(h.focus.y-h.origo.y)./(h.focus.x-h.origo.x).*(origin(:,1)-h.origo.x)+h.focus.y-h.focus.z*tan(h.tilt(2));
            end
            origin(:,3)=0;
            origin(isnan(origin))=0; % solve divisions by 0
            
            value=origin;
        end
        
        function value=get.propagation_distance(h)
            value=sqrt(sum(([h.focus.x-h.origin(:,1) h.focus.y-h.origin(:,2) h.focus.z-h.origin(:,3)]).^2,2));
        end
        
        function [ratio_theta ratio_phi]=incidence_ratios(h)
            if isempty(h.sequence) % no sequence then it must be a receive aperture
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                x=ones(h.focus.N_pixels,1)*(h.probe.x.'); 
                y=ones(h.focus.N_pixels,1)*(h.probe.y.'); 
                z=ones(h.focus.N_pixels,1)*(h.probe.z.'); 

                if isa(h.focus,'uff.sector_scan')
                    % distances
                    x_dist=h.focus.apex.x-x;
                    y_dist=h.focus.apex.y-y;
                    z_dist=h.focus.z*ones(1,h.probe.N_elements)-z;
                else
                    % distances
                    x_dist=h.focus.x*ones(1,h.probe.N_elements)-x;
                    y_dist=h.focus.y*ones(1,h.probe.N_elements)-y;
                    z_dist=h.focus.z*ones(1,h.probe.N_elements)-z;
                end
                
                % clamping z=0
                z_dist(abs(z_dist)<1e-6)=1e-6;

                % azimuth and elevation tangents
                tan_theta = x_dist./z_dist;
                tan_phi = y_dist./z_dist;

            else % transmit apodization
                for n=1:length(h.sequence)
                    % plane wave
                    if (h.sequence(n).wavefront==uff.wavefront.plane||isinf(h.sequence(n).source.distance))
                        z_dist(:,n)=h.focus.z;
                        % clamping z=0
                        z_dist(abs(z_dist)<1e-6)=1e-6;
                        tan_theta(:,n)=ones(h.focus.N_pixels,1)*h.sequence(n).source.azimuth;
                        tan_phi(:,n)=ones(h.focus.N_pixels,1)*h.sequence(n).source.elevation;
                    % diverging or converging waves
                    else
                        % distances
                        x_dist=h.focus.x-h.sequence(n).source.x;
                        y_dist=h.focus.y-h.sequence(n).source.y;
                        z_dist(:,n)=h.focus.z-h.sequence(n).source.z;
                        % clamping z=0
                        z_dist(abs(z_dist)<1e-6)=1e-6;
                        % azimuth and elevation tangents
                        tan_theta(:,n) = x_dist./z_dist(:,n);
                        tan_phi(:,n) = y_dist./z_dist(:,n);
                    end
                end
            end
            
            % ratios
            ratio_theta = abs(h.f_number(1)*tan_theta);
            ratio_phi = abs(h.f_number(2)*tan_phi);
            
            % minimum aperture
            min_theta_ratio=abs(z_dist.*tan_theta/h.minimum_aperture(1)); 
            min_phi_ratio=abs(z_dist.*tan_phi/h.minimum_aperture(2));
            min_theta_mask=ratio_theta>min_theta_ratio;
            min_phi_mask=ratio_phi>min_phi_ratio;
            ratio_theta(min_theta_mask)=min_theta_ratio(min_theta_mask);
            ratio_phi(min_phi_mask)=min_phi_ratio(min_phi_mask);
            
            % maximum aperture
            max_theta_ratio=abs(z_dist.*tan_theta/h.maximum_aperture(1)); 
            max_phi_ratio=abs(z_dist.*tan_phi/h.maximum_aperture(2));
            max_theta_mask=ratio_theta<max_theta_ratio;
            max_phi_mask=ratio_phi<max_phi_ratio;
            ratio_theta(max_theta_mask)=max_theta_ratio(max_theta_mask);
            ratio_phi(max_phi_mask)=max_phi_ratio(max_phi_mask);
        end
        
        function value=get.N_elements(h)
            if isempty(h.sequence)
                assert(numel(h.probe)>0,'The PROBE parameter is not set.');
                value=h.probe.N_elements; 
            else
                value=length(h.sequence);
            end
        
        end
        
        %% display methods 
        function figure_handle=plot(h,figure_handle_in,n)
            % PLOT Plots channel data
            if nargin>1 && not(isempty(figure_handle_in))
                figure_handle=figure(figure_handle_in);
            else
                figure_handle=figure();
            end
            
            if nargin <3
                n=round(size(h.data,2)/2);
            end
            
            colorMap = tools.inferno;
            
            isreceive = isempty(h.sequence);
            
            if isa(h.focus,'uff.linear_scan')
                imagesc(h.focus.x_axis*1e3,h.focus.z_axis*1e3,reshape(h.data(:,n),[h.focus.N_z_axis h.focus.N_x_axis]))
                xlabel('x [mm]');
                ylabel('z [mm]');
                set(gca,'Ydir','reverse');
                set(gca,'fontsize',14);
                colorbar;
                %caxis([0 1]);
                if isreceive
                    title(sprintf('Apodization values for element %d',n));
                else
                    title(sprintf('Apodization values for wave %d',n));
                end
            
            elseif isa(h.focus,'uff.sector_scan')
                 x_matrix=reshape(h.focus.x,[h.focus(1).N_depth_axis h.focus(1).N_azimuth_axis]);
                 z_matrix=reshape(h.focus.z,[h.focus(1).N_depth_axis h.focus(1).N_azimuth_axis]);
                 
                 
                 subplot(1,2,1);
                 pcolor(x_matrix*1e3,z_matrix*1e3,reshape(h.data(:,n),[h.focus.N_depth_axis h.focus.N_azimuth_axis]));
                 xlabel('x [mm]');
                 ylabel('z [mm]');
                 shading('flat');
                 set(gca,'fontsize',14);
                 set(gca,'YDir','reverse');
                 axis('tight','equal');
                 colorbar();
                 %caxis([0 1]);
                 if isreceive
                    title(sprintf('Apodization for element %d. Click or Enter.',n));
                 else
                    title(sprintf('Apodization for wave %d. Click or Enter.',n));
                 end
                 [x z]=ginput(1);
                 while ~isempty(x)
                     [~, ns]=min(sum((h.focus(1).xyz-[x 0 z]/1e3).^2,2));
                     subplot(1,2,2);
                     plot(h.data(ns,:)); grid on; axis tight;
                     ylim([0 1.2]);
                     if isreceive
                        title(sprintf('Receive apodization at pixel (%0.2f,%0.2f) mm.',x,z));
                        xlabel('Element');
                     else
                        title(sprintf('Transmit apodization at pixel (%0.2f,%0.2f) mm.',x,z));
                        xlabel('wave');
                     end
                     set(gca,'fontsize',14);
                     
                     subplot(1,2,1);
                     [x z]=ginput(1);
                 end
            else
                error('Only apodization plot for uff.linear_scan and uff.sector_scan are supported for now.');
            end
        end
    end
end