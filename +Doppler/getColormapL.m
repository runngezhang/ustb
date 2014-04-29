function map1 = getColormapL()
% colormap

map=zeros(256*2,3);

%positive red values

map(1:15,1)=linspace(1/255,128/255,15)';

map(16:128,1)=linspace(128/255,1,113)';

map(129:256,1)=ones(256-129+1,1);

%positive green values

map(1:15,2)=zeros(15,1);

map(16:128,2)=linspace(0,64/255,113)';

map(129:256,2)=linspace(64/255,1,256-129+1)';

%positive blue values

map(1:256,3)=zeros(256,1);

%negative red values

map(257:512,1)=zeros(256,1);

%negative green values

map(1+256:15+256,2)=zeros(15,1);

map(16+256:128+256,2)=linspace(0,64/255,113)';

map(129+256:256+256,2)=linspace(64/255,1,256-129+1)';

%negative blue values

map(1+256:15+256,3)=linspace(1/255,128/255,15)';

map(16+256:128+256,3)=linspace(128/255,1,113)';

map(129+256:256+256,3)=ones(256-129+1,1);

map(1:256,:) = flipud(map(1:256,:));

map1(:,1)=map(1:2:end,1);

map1(:,2)=map(1:2:end,2);

map1(:,3)=map(1:2:end,3);

map1=flipud(map1);