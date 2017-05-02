classdef das_cc_processes < beamformer
%DAS_CC_PROCESSES   Example of DAS and coherent compounding beamformer 
% using USTB built-in processes
%
%   See also BEAMFORMER

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Date: 2017/05/02$

    %% constructor
    methods (Access = public)
        function h=das_cc_processes(bmf)
            %DAS_CC_PROCESSES   Constructor of the das_cc_processes class
            %
            %   Syntax:
            %   h = das_cc_processes()
            %
            %   See also BEAMFORMER                      
           
            if nargin<1
                bmf=[]; 
            end
            
            h=h@beamformer(bmf);
            
            % Logistics
            h.name='DAS Cohe Comp processes';   
            h.reference= 'www.ustb.no';                
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};    
            h.version='v1.0.1';     
        end
    end
    
    
    %% go method
    methods 
        function out_data = go(h)
            % declare processes
            proc_1=process.das_matlab();
            proc_2=process.coherent_compounding();
               
            % input data to first process
            proc_1.channel_data=h.channel_data;
            proc_1.receive_apodization=h.receive_apodization;
            proc_1.transmit_apodization=h.transmit_apodization;
            proc_1.scan=h.scan;
                             
            % go process 1 & set input data to process 2
            proc_2.beamformed_data = proc_1.go();
                            
            % go process 2
            out_data = proc_2.go();                   
        end
    end    
end