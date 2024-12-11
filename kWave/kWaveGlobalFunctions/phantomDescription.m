function phantomConfirm = phantomDescription(codename)
if(strcmp(codename, 'CystSingleAnechoic'))
    disp('Selected Phantom: Single anechoic cyst phantom for contrast experiments');
end

if(strcmp(codename,'CystSingleHyperEchoic'))
    disp('Selected Phantom: Single hyperechoic cyst phantom for contrast experiments');
end

if(strcmp(codename,'PinMulti'))
    disp('Selected Phantom: Multiple pin for resolution experiments');
end

if(strcmp(codename,'specularReflLine'))
    disp('Selected Phantom: Specular reflector line');
end

if(strcmp(codename,'specularReflSupraSN'))
    disp('Selected Phantom: Musculoskeletal/ Normal supraspinatous tendon');
end

if(strcmp(codename,'SpeckleBackground'))
    disp('Selected Phantom: Speckle background');
end

if(strcmp(codename, 'LiverNormal'))
    disp('Selected Phantom: Normal liver phantom');
    %disp('Preferred frequency range: 2-5MHz');
end

if(strcmp(codename, 'LiverFatty'))
    disp('Selected Phantom: Fatty liver phantom');
    %disp('Preferred frequency range: 2-5MHz');
end

if(strcmp(codename, 'LiverFibrotic'))
    disp('Selected Phantom: Fibrotic liver phantom');
    %disp('Preferred frequency range: 2-5MHz');
end

if(strcmp(codename, 'LungNormalAlines') || strcmp(codename,'LungNormalAlinesImg'))
    disp('Selected Phantom: Normal lung phantom');
end

if(strcmp(codename, 'LungAbnormalBlines') || strcmp(codename,'LungAbnormalBlinesImg'))
    disp('Selected Phantom: Abnormal lung phantom -Blines');
end

if(strcmp(codename, 'LungAbnormalConsolidation') || strcmp(codename,'LungAbnormalConsolidationImg'))
    disp('Selected Phantom: Abnormal lung phantom - Plueural consolidation');
end



prompt = 'Phantom selection confirmed? Y/N [Y]:';
str = input(prompt,'s');
if ((str == 'Y') || (str == 'y'))
    phantomConfirm = 1;
elseif isempty(str)
    str = 'Y';
    phantomConfirm = 1;
else
    phantomConfirm = 0;
end
end