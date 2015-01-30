classdef dataset < handle
%DATASET    Ultrasound dataset.
%
%   DATASET is a superclass containing the properties and methods that are
%   common to a number of ultrasound datasets: Synthetic transmitaperture (sta), 
%   Coherent plane wave (cpw), Virtual source (vs). 
%
%   See also STA, CPW, VS

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    properties  (SetAccess = public)
        name=''             % String containing the name of the dataset
        creation_date       % String containing the date the dataset class was created
        format = E.signal_format.IQ  % Signal_format RF or IQ (enumerations.signal_format.RF, default=enumerations.signal_format.IQ)
        geom                % matrix M x 3 containing probe geometry [x, y, z] (m)
        data                % Matrix containing the numerical data. For acquisition
                            % datasets the matrix dimensions are [time_samples, channels, firings, frames]
                            % For image reconstruction datasets the matrix
                            % dimensions are [pixels, frames]
        time                % vector containing fast time (s)
        c0                  % reference speed of sound (m/s)
        modulation_frequency % value conatining the modulation frequency (Hz), only required for IQ format
    end
    
    properties  (SetAccess = protected)   
        frames                  % number of frames in the dataset        
        channels                % number of channels in the transducer
        firings                 % number of firings in the sequence (i.e. number of plane waves in CPWI, virtual sources in VSI, or elements in STAI)
        initial_time            % initial time (s)
        sampling_frequency      % sampling frequency (Hz)
        transmit_apodization    % matrix containing the apodization used for transmiting
        receive_apodization     % matrix containing the apodization used for receiving
    end
    
    %% constructor
    methods (Access = public)
        function h = dataset(name)
            %DATASET    Constructor of the DATASET class.
            %
            %   Syntax:
            %   DATASET(name) 
            %       name                    Name of the dataset
            %
            %   See also STA, CPW, DW
 
            if exist('name') h.name=name; end
            t_now=now;
            h.creation_date=[date sprintf('-%d-%d-%d',hour(t_now),minute(t_now),round(second(t_now)))];
        end
    end
    
    %% Relaunching method
    methods (Access = public)
        function sig=launch_implementation(h,r,imp)
            % Launches appropriate implementation
            switch(imp)
                case E.implementation.simple_matlab
                    sig=h.ir_simple_matlab(r);
                case E.implementation.optimized_matlab
                    sig=h.ir_optimized_matlab(r);
                case E.implementation.mex
                    sig=h.ir_mex(r);
                case E.implementation.mex_gpu
                    sig=h.ir_mex_gpu(r);
                otherwise
                    error('Selected implementation is not supported');
            end
        end
    end
    
    %% implementation methods, to be overloaded
    methods 
        function sig=ir_simple_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Simple matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_optimized_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Optimazed matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex_gpu(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex gpu: The implementation is not available yet!');
            sig=[];
        end
    end
    
    %% set methods
    methods  
        function set.time(h,input_time)
            assert(size(input_time,1)>size(input_time,2), 'The time vector must be a column vector!')
            h.time=input_time;
            h.sampling_frequency=1/mean(diff(input_time));
            h.initial_time=input_time(1);
        end
        function set.data(h,input_data)
            assert(ndims(input_data)>=3, 'The dataset must have 4 dimensions [time_samples, channels, firings, frames]')
            h.channels=size(input_data,2);
            h.firings=size(input_data,3);
            h.frames=size(input_data,4);
            h.data=input_data;
        end
        function set.geom(h,input_geom)
            assert(ndims(input_geom)==2&&size(input_geom,2)==3, 'Wrong probe geometry definition. It should be a three column matrix [x, y, z]')
            h.geom=input_geom;
        end
    end
    
end

