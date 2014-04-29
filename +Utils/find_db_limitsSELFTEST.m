function find_db_limitsSELFTEST()

y = [0:0.25:1,0.75:-0.25:0];
y = [y,fliplr(y)];
x = 1:length(y);
lim = -6;
[db_ix db_x] = find_db_limits(y,lim,x);
assert(db_ix(1) == 4 && db_ix(2) == 16,'HalfWave:selfTestFailed','find_db_limits found wrong limits');

% figure(1)
% clf
% plot(x,y,db_x([1 1]),[0 1],'k--',db_x([2 2]),[0 1],'k--');

tight = 1;
[db_ix db_x] = find_db_limits(y,lim,x,tight);
assert(db_ix(1) == 4 && db_ix(2) == 8,'HalfWave:selfTestFailed','find_db_limits found wrong limits');

% hold on
% plot(x,y,db_x([1 1]),[0 1],'r-.',db_x([2 2]),[0 1],'r-.');
% hold off

y2 = y;
y2(5:end) = 1;
[db_ix db_x] = find_db_limits(y2,lim,x,0);
assert(db_ix(1) == 4 && db_ix(2) == 18,'HalfWave:selfTestFailed','find_db_limits found wrong limits');

y3 = y;
y3(:) = 1;
[db_ix db_x] = find_db_limits(y3,lim,x,0);
assert(db_ix(1) == 1 && db_ix(2) == 18,'HalfWave:selfTestFailed','find_db_limits found wrong limits');
