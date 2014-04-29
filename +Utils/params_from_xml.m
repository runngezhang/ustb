function params_struct = params_from_xml(filename)

params_struct = [];
try 
    xmlDoc = xmlread(filename);
catch
    disp('could not read xml file');
    params_struct = [];
    return;
end;

xmlUlts = xmlDoc.getElementsByTagName('parameters');
if xmlUlts.getLength() == 0
    params_struct = [];
    return;
end

xmlUlt = xmlUlts.item(0);
xmlParameters = xmlUlt.getElementsByTagName('parameter');

for k=0:xmlParameters.getLength()-1
    paramNode = xmlParameters.item(k);        
    valueNodes = paramNode.getElementsByTagName('value');
    
    id = str2double(char(paramNode.getAttribute('id')));
    typeid = str2double(char(paramNode.getAttribute('typeid')));
    name = char(paramNode.getAttribute('name'));
    unit = char(paramNode.getAttribute('unit'));
    
    params_struct{id}.id = id;
    params_struct{id}.name = name;
    params_struct{id}.typeid = typeid;
    params_struct{id}.unit = unit;
    
    if valueNodes.getLength > 0
        valueNode = valueNodes.item(0);
        while ~isempty(valueNode)
            if valueNode.getNodeType() == valueNode.ELEMENT_NODE                                
                switch typeid
                    case 0
                        params_struct{id}.value = str2double(char(valueNode.getFirstChild.getData()));                        
                    case 2     
                        nodes = valueNode.getChildNodes();                        
                        params_struct{id}.value = struct('left',0,'right',0,'bottom',0,'top',0);
                        for ii=0:(nodes.getLength()-1)
                            node = nodes.item(ii);
                            if node.getNodeType() == node.ELEMENT_NODE
                                value = str2double(char(node.getFirstChild.getData()));
                                nodeName = char(node.getNodeName());                                
                                params_struct{id}.value.(nodeName) = value;                                                                
                            end
                        end                        
                    case 6
                        nodes = valueNode.getChildNodes();                        
                        params_struct{id}.value = struct('t',0,'m',0,'b',0,'vm',0);
                        for ii=0:(nodes.getLength()-1)
                            node = nodes.item(ii);
                            if node.getNodeType() == node.ELEMENT_NODE
                                value = str2double(char(node.getFirstChild.getData()));
                                nodeName = char(node.getNodeName());                                
                                params_struct{id}.value.(nodeName) = value;                                                                
                            end
                        end                        
                    case 8
                        params_struct{id}.value = zeros(1,8);
                        nodes = valueNode.getChildNodes();
                        for ii=0:(nodes.getLength()-1)
                            node = nodes.item(ii);
                            if node.getNodeType() == node.ELEMENT_NODE
                                nodeName = char(node.getNodeName());
                                indx = str2double(nodeName(end));
                                if ~isnan(indx)
                                    params_struct{id}.value(indx) = str2double(char(node.getFirstChild.getData()));
                                end
                            end
                        end
                end
            end
            valueNode = valueNode.getNextSibling;        
        end
    end
end