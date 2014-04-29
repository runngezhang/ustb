function success = load_scanner_setup(xml_fname,host)
% success = load_scanner_setup(xml_fname,host)
%
% Sets up the scanner according to the setup save in the xml file. If no
% xml file is sepcified a gui file chooser is opened and you can choose
% your file.

success = 0;
if nargin < 1 || isempty(xml_fname)
    [xml_fname,path_name] = uigetfile('*.xml','Choose xml file');
    if xml_fname == 0
        return;
    end
    xml_fname = fullfile(path_name,xml_fname);
end

if nargin < 2
    host = 'localhost';
end

params = params_from_xml(xml_fname);
success = 1;
if ~isempty(params)
    ult = multerius;
    if ~ult.connect(host)
        return;
    end
    
    drop_params = [22 23 26 27 103 104 822 1121 1125 1133 1156 1157 1207 1361];
    for k=1:length(params)
        if isempty(params{k})
            continue;
        end
        
        if any(drop_params == params{k}.id)
            continue;
        end
        disp(params{k}.name);        
        switch params{k}.typeid
            case 0
                if ~ult.set_int_param(params{k}.id,params{k}.value)
                    warning('Could not set parameter %s with ID = %d',params{k}.name,params{k}.id);
                    success = 0;
                end
            case 2               
%                 if ~ult.set_rect_param(params{k}.id,params{k}.value)
%                     warning('Could not set parameter %s with ID = %d',params{k}.name,params{k}.id);
%                     success = 0;
%                 end
            case 6
                if ~ult.set_curve_param(params{k}.id,params{k}.value)
                    warning('Could not set parameter %s with ID = %d',params{k}.name,params{k}.id);
                    success = 0;
                end
        end               
    end
    ult.disconnect();    
end