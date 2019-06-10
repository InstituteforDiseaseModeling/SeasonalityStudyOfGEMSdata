clear all;
close all;

data = readtable('IntermediateProcessingMatlabDataStructure/estimated-rotavirus-fortnightly.csv');

%% Putting Data into format for analysis

% Put date string into datenum
dates = datenum(data.date);

% Estimated Case load
cases = zeros(length(dates),1);

for dateIndex = 1 : length(dates)
    
    % Seems to be some NA entries in the file - Need to talk with Dennis
    if isempty(str2num(data.estimatedRotavirusPositiveEligible{dateIndex}))
        cases(dateIndex,1) = 0;
    else
        cases(dateIndex,1) = str2num(data.estimatedRotavirusPositiveEligible{dateIndex});
    end
end

% Place
uniqueCountryNames = unique(data.country);
countryVec = zeros(length(dates),1);

for countryIndex = 1 : length(uniqueCountryNames)
   countryVec = countryVec + countryIndex * strcmp(uniqueCountryNames(countryIndex),data.country); 
end

% Clear up some memory
clear data dateIndex countryIndex;

%% Compute PowerSpectrum for each country

for countryIndex = 1 : 1 
   
    tempCaseData = cases(countryVec == countryIndex);
    tempDateData = dates(countryVec == countryIndex);
       
    h3 = figure;
    
    Pfa = [50 10 1 0.1]/100;
    Pd = 1-Pfa;
    [pxx,f,pth] = plomb(tempCaseData,tempDateData/365,'normalized','Pd',Pd);

    plot(f,pxx,'k','LineWidth',1.5), hold on
    h1 = plot(f,pth*ones(size(f')));
    %legend(h1,'FA 50%','FA 10%','FA 1%','FA 0.1%')
    axis([f(1) 6 min(pxx) max(pxx)])
    set(gca,'FontSize',22,'LineWidth',1.5)
    set(gca,'XTick',1:6,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    %xlabel('Frequency ~ 1/Year')
    %title('Lomb-Scargle periodogram')
    %axis([0 0.12 0 35]) 
    axis square
    set(gcf,'Position',[24 5 500 500]);
    set(gcf,'PaperPositionMode','auto')
    %saveas(h3,['figures\Analysis',uniqueCountryNames{countryIndex}],'pdf');
    %print(h3,['figures\Analysis',uniqueCountryNames{countryIndex}]);
    %savefig(h3,['figures\Analysis',uniqueCountryNames{countryIndex}],'pdf');
    print('-depsc2', '-loose',['figures\PLOSNTDPaper\Analysis',uniqueCountryNames{countryIndex}]);
    %close(h3);

end