classdef beamformed_data < handle
    %BEAMFORMED_DATA   beamformed_data definition. Children of HANDLE class
    %
    %   See also PULSE, PHANTOM, PROBE
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal (olemarius@olemarius.net)
    %
    %   $Last updated: 2017/03/20 $
    
    %% compulsory properties
    properties  (SetAccess = public)
        scan                       % SCAN class
        data                       % data
    end
 
    %% Logistics 
    properties  (SetAccess = public)
        name={}              % name of the dataset
        reference={}         % reference to the publication where it was used/acquired
        author={}            % contact of the authors
        version={}           % version of the dataset
    end
    
    %% optional properties
    properties  (SetAccess = public)
        phantom                    % PHANTOM class [optional]
        wave                       % WAVE class [optional]
        probe                      % PROBE class [optional]
        pulse                      % PULSE class [optional]
        beamformer                 % String with the beamformer specification [optional]
       sampling_frequency         % Sampling frequency in depth / z-direction        
    end
    
    %% dependent properties
    properties  (Dependent)
        N_pixels                   % number of pixels  
     end
    
    %% private properties
    properties (Access = private)
       current_frame             % If multi frame, this is the current frame shown 
       all_images
       image_handle
       in_title
       play_loop
    end
    
    %% constructor
    methods (Access = public)
        function h=beamformed_data(in_beamformed_data)
            %BEAMFORMED_DATA   Constructor of beamformed_data class
            %
            %   Syntax:
            %   h = beamformed_data()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE
            if nargin>0 && ~isempty(in_beamformed_data)
                h.copy(in_beamformed_data);
            end
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
                    eval(sprintf('h.%s = object.%s;',property_name,property_name));
                end
            end
        end
    end
    
    %% plot methods
    methods (Access = public)
        function figure_handle=plot(h,figure_handle_in,in_title,dynamic_range,compression)
            %PLOT Plots beamformed data
            %
            % Usage: figure_handle=plot(figure_handle,title,dynamic_range)
            %
            %   figure_handle   Handle to the figure to plot to (default: none)
            %   title           Figure title (default: none)
            %   dynamic_range   Displayed dynamic range (default: 60 dB)
            %   compression     String specifying compression type: 'log','none','sqrt' (default: 'log')
            
            if (nargin>1 && ~isempty(figure_handle_in) && isa(figure_handle_in,'matlab.ui.Figure')) || ...
                    (nargin>1 && ~isempty(figure_handle_in) && isa(figure_handle_in,'double'))
                figure_handle=figure(figure_handle_in);
                axis_handle = gca(figure_handle);
                hold on;
            elseif nargin>1 && ~isempty(figure_handle_in) && isa(figure_handle_in,'matlab.graphics.axis.Axes')
                axis_handle = figure_handle_in;
            else
                figure_handle=figure();
                axis_handle = gca(figure_handle);
            end

            if nargin<3
                h.in_title='Beamformed data';
            else
                h.in_title = in_title;
            end
            if nargin<4
                dynamic_range=60;
            end
            if nargin<5
                compression='log';
            end
            %Draw the image
            h.draw_image(axis_handle,h.in_title,dynamic_range,compression);
            
            %If more than one frame, add the GUI buttons
            if size(h.data,2) > 1 
                set(figure_handle, 'Position', [100, 100, 600, 700]);
                h.current_frame = 1;
                h.add_buttons(figure_handle);
                h.play_loop = 0;
                title([h.in_title,', Frame = ',num2str(h.current_frame),'/',num2str(size(h.all_images,3))]);
            end
        end
        
        function draw_image(h,axis_handle,in_title,dynamic_range,compression)
            
            % compress values
            switch compression
                case 'log'
                    envelope=abs(h.data);
                    envelope=20*log10(envelope./max(envelope(:)));
                    max_value=0;
                    min_value=-dynamic_range;
                case 'sqrt'
                    envelope=sqrt(abs(h.data));
                    max_value=max(envelope(:));
                    min_value=10^(-dynamic_range/20);
                case 'none'
                    envelope=abs(h.data);
                    max_value=max(envelope(:));
                    min_value=10^(-dynamic_range/20);
            end
                
            switch class(h.scan)
                case 'uff.linear_scan'
                    x_matrix=reshape(h.scan.x,[h.scan.N_z_axis h.scan.N_x_axis]);
                    z_matrix=reshape(h.scan.z,[h.scan.N_z_axis h.scan.N_x_axis ]);
                    h.all_images = reshape(envelope,[h.scan.N_z_axis h.scan.N_x_axis size(h.data,2)]);
                    h.image_handle = pcolor(axis_handle,x_matrix*1e3,z_matrix*1e3,h.all_images(:,:,1));
                    shading(axis_handle,'flat');
                    set(axis_handle,'fontsize',14);
                    set(axis_handle,'YDir','reverse');
                    axis(axis_handle,'tight','equal');
                    colorbar(axis_handle);
                    colormap(axis_handle,'gray');
                    xlabel(axis_handle,'x[mm]'); ylabel(axis_handle,'z[mm]');
                    caxis(axis_handle,[min_value max_value]);
                    title(axis_handle,in_title);
                    drawnow;
                case 'uff.linear_3D_scan'
                    [radial_matrix axial_matrix] = meshgrid(h.scan.radial_axis,h.scan.axial_axis);
                    h.all_images = reshape(envelope,[h.scan.N_axial_axis h.scan.N_radial_axis size(h.data,2)]);
                    [az,el] = view();
                    if (el==90) 
                        % plot in 2D
                        h.image_handle = pcolor(axis_handle,radial_matrix*1e3,axial_matrix*1e3,h.all_images(:,:,1));
                        shading(axis_handle,'flat');
                        set(axis_handle,'fontsize',14);
                        set(axis_handle,'YDir','reverse');
                        axis(axis_handle,'tight','equal');
                        colorbar(axis_handle);
                        colormap(axis_handle,'gray');
                        xlabel(axis_handle,'radial[mm]'); ylabel(axis_handle,'axial[mm]');
                        caxis(axis_handle,[min_value max_value]);
                        title(axis_handle,in_title);
                    else
                        % plot in 3D
                        x_matrix=reshape(h.scan.x,[h.scan.N_axial_axis h.scan.N_radial_axis]);
                        y_matrix=reshape(h.scan.y,[h.scan.N_axial_axis h.scan.N_radial_axis]);
                        z_matrix=reshape(h.scan.z,[h.scan.N_axial_axis h.scan.N_radial_axis]);
                        surface(axis_handle);
                        surface(x_matrix*1e3,y_matrix*1e3,z_matrix*1e3,h.all_images(:,:,1));
                        shading(axis_handle,'flat');
                        set(axis_handle,'fontsize',14);
                        %set(axis_handle,'YDir','reverse');
                        axis(axis_handle,'tight','equal');
                        colorbar(axis_handle);
                        colormap(axis_handle,'gray');
                        xlabel(axis_handle,'x[mm]'); 
                        ylabel(axis_handle,'y[mm]');
                        zlabel(axis_handle,'z[mm]');
                        caxis(axis_handle,[min_value max_value]);
                        title(axis_handle,in_title);                        
                    end
                    drawnow;
                case 'uff.sector_scan'
                    x_matrix=reshape(h.scan.x,[h.scan.N_depth_axis h.scan.N_azimuth_axis]);
                    z_matrix=reshape(h.scan.z,[h.scan.N_depth_axis h.scan.N_azimuth_axis ]);
                    h.all_images = reshape(envelope,[h.scan.N_depth_axis h.scan.N_azimuth_axis size(h.data,2)]);
                    h.image_handle = pcolor(axis_handle,x_matrix*1e3,z_matrix*1e3,h.all_images(:,:,1));
                    shading(axis_handle,'flat');
                    set(axis_handle,'fontsize',14);
                    set(axis_handle,'YDir','reverse');
                    axis(axis_handle,'tight','equal');
                    colorbar(axis_handle);
                    colormap(axis_handle,'gray');
                    xlabel(axis_handle,'x[mm]'); ylabel(axis_handle,'z[mm]');
                    caxis(axis_handle,[min_value max_value]);
                    title(axis_handle,in_title);
                    drawnow;
                otherwise
                    error(sprintf('Dont know how to plot on a %s yet. Sorry!',class(b_data.scan)));
            end
        end
        
        function img = get_image(h)
            envelope=abs(h.data);
            envelope=20*log10(envelope./max(envelope(:)));
            switch class(h.scan)
                case 'uff.linear_scan'
                    img = reshape(envelope,[h.scan.N_z_axis h.scan.N_x_axis size(h.data,2)]);
                case 'uff.sector_scan'
                    img = reshape(envelope,[h.scan.N_depth_axis h.scan.N_azimuth_axis size(h.data,3)]);
                otherwise
                    error(sprintf('Dont know how to plot on a %s yet. Sorry!',class(b_data.scan)));
            end
        end
        
        function h = calculate_sampling_frequency(h,c)
            assert(exist('c')==1,'Please give speed of sound as input');
            % calculate sampling frequency
            if isa(h.scan,'uff.linear_scan')
               dz = h.scan.z_step;
            elseif isa(h.scan,'uff.sector_scan')
               dz = h.scan.depth_step;
            end  
            h.sampling_frequency = (c/dz/2); % effective sampling frequency (Hz)      
        end
    end
    
    %% set methods
    methods
        function h=set.phantom(h,in_phantom)
            if ~isempty(in_phantom)
                assert(isa(in_phantom,'uff.phantom'), 'The input is not a PHANTOM class. Check HELP PHANTOM.');
                h.phantom=in_phantom;
            end
        end
        function h=set.pulse(h,in_pulse)
            if ~isempty(in_pulse)
                assert(isa(in_pulse,'uff.pulse'), 'The input is not a PULSE class. Check HELP PULSE.');
                h.pulse=in_pulse;
            end
        end
        function h=set.probe(h,in_probe)
            if ~isempty(in_probe)
                assert(isa(in_probe,'uff.probe'), 'The input is not a PROBE class. Check HELP PROBE.');
                h.probe=in_probe;
            end
        end
        function h=set.wave(h,in_wave)
            if ~isempty(in_wave)
                assert(isa(in_wave,'uff.wave'), 'The input is not a WAVE class. Check HELP WAVE.');
                h.wave=in_wave;
            end
        end
        function h=set.scan(h,in_scan)
            if ~isempty(in_scan)
                assert(isa(in_scan,'uff.scan'), 'The input is not a SCAN class. Check HELP SCAN.');
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
        function value=get.N_pixels(h)
            value=h.scan.N_pixels;
        end        
    end
    
    %% GUI functions
    methods (Access = private)
        function add_buttons(h,figure_handle)
            ButtonH=uicontrol('Parent',figure_handle,'Style','pushbutton','String','Previous frame','Units','normalized','Position',[0.22 0.95 0.2 0.05],'Visible','on','Callback',{@h.plot_previous_frame,h});
            ButtonH=uicontrol('Parent',figure_handle,'Style','pushbutton','String','Play movie loop','Units','normalized','Position',[0.42 0.95 0.2 0.05],'Visible','on','Callback',{@h.play_movie_loop,h});
            ButtonH=uicontrol('Parent',figure_handle,'Style','pushbutton','String','Next frame','Units','normalized','Position',[0.62 0.95 0.2 0.05],'Visible','on','Callback',{@h.plot_next_frame,h});
        end
        
        function plot_previous_frame(h,var1,var2,var3)
            if h.current_frame-1 > 0
                h.current_frame = h.current_frame - 1;
                set(h.image_handle,'CData',h.all_images(:,:,h.current_frame));
                title([h.in_title,', Frame = ',num2str(h.current_frame),'/',num2str(size(h.all_images,3))]);
            end
        end
        
        function plot_next_frame(h,var1,var2,var3)
            if h.current_frame+1 <= size(h.all_images,3)
                h.current_frame = h.current_frame + 1;
                set(h.image_handle,'CData',h.all_images(:,:,h.current_frame));
                title([h.in_title,', Frame = ',num2str(h.current_frame),'/',num2str(size(h.all_images,3))]);
            end
        end
        
        function play_movie_loop(h,var1,var2,var3)
            h.play_loop = ~h.play_loop;
            while(h.play_loop)
                h.current_frame = h.current_frame+1;
                if h.current_frame > size(h.all_images,3)
                    h.current_frame = 1;
                end
                set(h.image_handle,'CData',h.all_images(:,:,h.current_frame));
                title([h.in_title,', Frame = ',num2str(h.current_frame),'/',num2str(size(h.all_images,3))]);
                drawnow();
                pause(0.05);
            end
        end
    end
    
    % get object name
    methods
    	function  out  = objname(h)
            out = evalin('caller','inputname(1)');
        end
    end
    
    %% Ultrasound File Format (UFF)
    methods (Access = public)    
        function write(h, uff_object, location)
            assert(isa(uff_object,'uff'),'The input is not a UFF object')
            if nargin<3 location=[]; end

            fprintf('UFF write: %s -> %s',h.objname,uff_object.filename);
            tic;
            uff_object.write(location,h,h.objname);
            fprintf(' [%0.1fs]\n',toc);
            
        end
    end
end