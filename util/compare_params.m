function diff_params = compare_params(params1,params2)
cnt = 1;
for kk=1:length(params1)    
    if isempty(params1{kk})
        continue;
    end    
    if params1{kk}.typeid == 0
        if params1{kk}.value ~= params2{kk}.value            
            fprintf('Param Name (ID) : %s (%d)\n',params1{kk}.name,params1{kk}.id);
            fprintf('Param value: %d :: %d\n',params1{kk}.value,params2{kk}.value);
            disp('---------------------------')
            diff_params{cnt}.p1 =  params1{kk};
            diff_params{cnt}.p2 =  params2{kk};
            cnt = cnt + 1;
        end
    elseif params1{kk}.typeid == 6
        v1 = [params1{kk}.value.t,params1{kk}.value.m,params1{kk}.value.b,params1{kk}.value.vm];
        v2 = [params2{kk}.value.t,params2{kk}.value.m,params2{kk}.value.b,params2{kk}.value.vm];
        if v1(1) ~= v2(1) || v1(2) ~= v2(2) || v1(3) ~= v2(3) || v1(4) ~= v2(4)            
            fprintf('Param Name (ID) : %s (%d)\n',params1{kk}.name,params1{kk}.id);
            fprintf('Param value: [%d,%d,%d,%d] :: [%d,%d,%d,%d]\n',v1(1),v1(2),v1(3),v1(4),v2(1),v2(2),v2(3),v2(4));
            disp('---------------------------')
            diff_params{cnt}.p1 =  params1{kk};
            diff_params{cnt}.p2 =  params2{kk};
            cnt = cnt + 1;
        end
    end
    
end