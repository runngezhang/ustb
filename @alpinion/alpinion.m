classdef alpinion < handle
    %ALPINION   Class for reading Alpinion research scanner data to USTB
    
    %   authors: Ole Marius Hoel Rindal (olemarius@olemarius.net)
    %   $Date: 2017/02/17$
    
    properties (SetAccess = public)
        parameter_filename                    %Name on the file with parameters to read
        data_filename                         %Name on the file with data to read
    end
    
    %% Constructor
    methods (Access = public)
        function h = alpinion(parameter_filename,data_filename)
            %ALPINION Constructor of the Alpinion class
            %
            if exist('parameter_filename')h.parameter_filename = parameter_filename;end
            if exist('data_filename')h.data_filename = data_filename;end
        end
    end
    
    %% Read parameters
    methods (Access = public)
        function dataset = read(h,dataset)
            
            if ismember(class(dataset),'cpw')
                % Load the filename
                load(h.parameter_filename);
                load(h.data_filename);
                
                % Reading paramteres
                dataset.name = ['Alpinion Research Scanner, CPW, filename = ',h.data_filename];
                %dataset.sampling_frequency = 
                fs = double(System.Parameters.sampleFreqMHz*10^6); % sampling frequency in Hz
                dataset.c0 = double(Parameter{1}.speedOfSoundMps);
                dataset.angle = deg2rad(double(Roi{1}.steerAngleDegA));
                dataset.format = E.signal_format.RF;
                dataset.center_frequency = double(System.Transducer.frequencyMHz*10^6); %center frequency in Hz 
                
                % Calculate element positions
                nbr_elements = double(System.Transducer.elementCnt);
                element_pitch_mm = double(System.Transducer.elementPitchCm)/100;
                elpos = zeros(3,nbr_elements);
                elpos(1,:) = linspace(-element_pitch_mm*(nbr_elements/2),element_pitch_mm*(nbr_elements/2),nbr_elements);
                dataset.geom = elpos';
                
                % Read data
                var_names = who('AdcData_frame*'); % Names on variables containing RF data
                data_initial = eval(var_names{1});
                dt = 1/fs;  
                dataset.time = [0:dt:(size(data_initial,2)-1)*dt]'; %Set og initial time        

                CPW = zeros(size(dataset.time,1),size(dataset.geom,1),size(dataset.angle,1),1);
                for transmission = 1:length(var_names)
                    %Read RF data
                    data = eval(var_names{transmission});
                    rfData = double([data(:,:,1)' data(:,:,2)'...
                        data(:,:,3)' data(:,:,4)']);
                   % for elmt = 1:size(rfData,2)
                   %     rfDataAnalytical(:,elmt) = hilbert(rfData(:,elmt));
                   % end
                   % rfDataAllFramesAnalytical(:,:,transmission) = rfDataAnalytical;
                    
                   %Compensate for steering regarding t0
                    D = abs(dataset.geom(1,1)-dataset.geom(end,1));
                    q = abs((D/2)*sin(dataset.angle(transmission)));
                    t_in=dataset.time-q/(dataset.c0);
                    v_aux=interp1(t_in,rfData,dataset.time,'linear',0);
                    
                    % build the dataset
                    CPW(:,:,transmission)=v_aux;
                end
                
                dataset.data = CPW;
                    
            
            else
                error('Only CPW dataset is implemented for Alpinion');
            end
        end
    end
        
    
end

