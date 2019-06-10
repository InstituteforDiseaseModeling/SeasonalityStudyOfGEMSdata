% disease x country study of frequencies
% min k roh

addpath('matlabSupportingFunctions')

clc
%% process data again
retry = 1;

%% save figures
saveFig = 0;

if retry
    %% periodogram-related variables
    % probability values to test for false positive
    percs = [10 5 1 0.01];
    Pfa = percs/100;
    Pd  = 1-Pfa;
    % signal threshold, 3 -> 1%
    sth_ind = 1;
    
    %% figure-related variables
    % font size
    fs           = 14;
    legend_label = {['FA ',num2str(percs(1)),' %'],['FA ',num2str(percs(2)),' %'],['FA ',num2str(percs(3)),' %'],['FA ',num2str(percs(4)),' %']};
    
    %% Putting Data into format for analysis
    % 'country'    'fortnight'    'date'    'eligible'
    % columns 5 to 18 contain diseases data
    data            = readtable('../Data/estimated-positive-fortnightly.csv');
    diseaseNames    = data.Properties.VariableNames(5:end);
    all_eligible    = data.Properties.VariableNames(4);
    len_allDiseases = length(diseaseNames);
    
    % Put date string into datenum
    dates = datenum(data.date);
    
    % Length of total data
    len_data = length(dates);
    
    % Country names and their corresponding entries on the table
    %     'Bangladesh'
    %     'India'
    %     'Kenya'diseaseNames
    %     'Mali'
    %     'Mozambique'
    %     'Pakistan'
    %     'The Gambia'
    
    fprintf('*** Processing data with %d%% significance level ***\n',100-percs(sth_ind));
    
    countryNames = unique(data.country);
    cnLen        = length(countryNames);
    countryVec   = zeros(len_data,1);
    for countryIndex = 1 : cnLen
        countryVec = countryVec + countryIndex * strcmp(countryNames(countryIndex),data.country);
    end
    
    
    % all eligible children by country
    all_eligibleXcountry     = zeros(cnLen,1);
    % Cumulative incidence per disease
    diseaseIncidence         = zeros(cnLen,len_allDiseases);
    orderedDiseaseIndices    = zeros(cnLen,len_allDiseases);
    orderedDiseaseIncidences = zeros(cnLen,len_allDiseases);
    for countryIndex = 1 : cnLen
        fprintf('Obtaining cumulative disease incidence for %s\n', countryNames{countryIndex});
        incidenceSum                             = sum(data{countryIndex==countryVec,5:end});
        diseaseIncidence(countryIndex,:)         = incidenceSum;
        [sortedIS, indices]                      = sort(incidenceSum, 'descend');
        orderedDiseaseIndices(countryIndex,:)    = indices;
        orderedDiseaseIncidences(countryIndex,:) = sortedIS;
        all_eligibleXcountry(countryIndex)       = sum(data{countryIndex==countryVec,4});
    end
    
    
    % picked indices, sorted by total incidence
    % cutoff -  diseases with incidence where each country has at least 25 data points for each disease 
    % ranked by total incidence summed over all countries and all dates 
    % 0 index denotes all eligible children
    % Rankings:
    %     1.	EAEC	               10.	Campylobacter Coli or Jejuni
    %     2.	Rotavirus	           11.	Adenovirus
    %     3.	Giardia	               12.	Entamoeba
    %     4.	ETEC ST                13.	Aeromonas
    %     5.	Norovirus GII          14.	V cholerae
    %     6.	Shigella              
    %     7.	Cryptosporidium    
    %     8.	tEPEC                  
    %     9.	Norovirus GI           



    % Order of diseases in data file
    % (1) Aeromonas  (2) Campylobacter Coli or Jejuni (3) Shigella (4) V cholerae (5) ETEC ST (6) tEPEC (7) EAEC (8) Cryptosporidium (9) Rotavirus (10) Adenovirus 40/41 (11) Norovirus GI (12) Norovirus GII
    
    chosenInds = [0, 2, 3, 5, 6, 8, 9, 10]+4;
    % sorted inds: 7  9  5  12  3  8  6  11  2  10  1  4
    chosenDN   = {'Eligible', 'C. jejuni-coli', 'Shigella', 'ETEC ST', 'tEPEC', 'Cryptosporidium',  'Rotavirus', 'Adenovirus 40-41'};
    len_dis    = length(chosenInds);
    
    
    %% phase information
    phaseInfo  = zeros(cnLen, len_dis);
    sigInfo    = zeros(cnLen, len_dis); % store statistical significance 
    periodInfo = zeros(cnLen, len_dis); % store cycle length 
    
    %% data structure to store pxx min&max, f min&max for each disease
    axes_data = zeros(len_dis, 4);
    % set min to large value
    axes_data(:,[1 3]) = intmax;
    
    %% get a periodogram for each country and store the information
    for countryIndex = 1:cnLen
        fprintf('\n%s\n',countryNames{countryIndex})
        
        % date index for a particular country
        date_inds        = find(countryVec == countryIndex);
        len_di           = length(date_inds);
        % storage variable
        casesData        = zeros(len_dis,len_di);
        % dates
        dateData         = dates(date_inds);
        % check for NA entries
        for dateIndex = 1 : len_di
            for diseaseIndex = 1:len_dis
                if isempty(data{date_inds(dateIndex),chosenInds(diseaseIndex)})
                    casesData(diseaseIndex,dateIndex) = 0;
                else
                    casesData(diseaseIndex,dateIndex) = data{date_inds(dateIndex),chosenInds(diseaseIndex)};
                end
            end
        end
        % error checking, each entry should be > 50, i.e., filtering criterion
        % sum(casesData,2)'
        
        %% structure to hold information
        eval([regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind),' = struct;']);
        
        %% plot disease data
        axis_lim = [min(dateData) max(dateData)];
        for diseaseIndex = 1:length(chosenDN)
            dis_name = chosenDN{diseaseIndex};
            fprintf('\t%s\n',dis_name);
            caseDatai = casesData(diseaseIndex,:)';
            
            h = figure(1);
            set(h, 'Position', [100 80 636 800])
            subplot(3,1,1), plot(dateData,caseDatai,'k','LineWidth',1.5)
            axis([axis_lim 0 max(1,max(caseDatai))])
            ylabel([dis_name, ' Cases'])
            title([countryNames{countryIndex},' ',dis_name])
            set(gca,'FontSize',fs,'LineWidth',1.5)
            datetick('x','yyyy','keeplimits')
            [f,mx1,~] = powerSpectrum(caseDatai);
            fscaled = f/((dateData(end)-dateData(1))/365);
            
            subplot(3,1,2), loglog(fscaled, mean(mx1,2), 'Color', 'k');
            axis([fscaled(1) fscaled(end) min(mx1) max(1,max(mx1))])
            xlabel('Frequency ~ 1/Year')
            title('Power Spectrum')
            set(gca,'FontSize',fs,'LineWidth',1.5)
            
            % significance testing
            [pxx,f,pth] = plomb(caseDatai,dateData/365,'normalized','Pd',Pd);
            subplot(3,1,3), plot(f,pxx,'-ok')
            hold on
            h1 = plot(f,pth*ones(size(f')));
            legend(h1, legend_label)
            axis([f(1) f(end) max(0,min(pxx)) max(1,max(pxx))])
            set(gca,'FontSize',fs,'LineWidth',1.5)
            xlabel('Frequency ~ 1/Year')
            title('Lomb-Scargle periodogram')
            hold off
            % find frequencies that are above 's_th' threshold
            th_inds = pxx > pth(sth_ind);
            th_pds  = pxx(th_inds);
            th_f    = f(th_inds);
            if isempty(th_pds)
                th_pds = {};
            end
            % check annual cycle's significance
            % pick the tallest magnitude within 11 and 13 months then compare it with Pd threshold values
            % 1 - not significant
            % 2 - 10%
            % 3 - 5%
            % 4 - 1%
            % 5 - 0.1%
            annual_peak_indices = find(f>11/12 & f<13/12); % between 11 and 13 months
            annual_peaks        = pxx(annual_peak_indices);
            [max_ap, m_ind]     = max(annual_peaks);
            periodInfo(countryIndex, diseaseIndex) = f(annual_peak_indices(m_ind));
            max_ind = find(max_ap > pth,1,'last');
            if isempty(max_ind)
                sigInfo(countryIndex,diseaseIndex) = 1;
            else
                sigInfo(countryIndex,diseaseIndex) = max_ind+1; % magnitude at the phase
            end
            
            %% update axes data
            pxx_min = axes_data(diseaseIndex,1);
            pxx_max = axes_data(diseaseIndex,2);
            f_min   = axes_data(diseaseIndex,3);
            f_max   = axes_data(diseaseIndex,4);
            
            if ~isempty(th_pds)
                if pxx_min > min(th_pds)
                    axes_data(diseaseIndex,1) = min(th_pds);
                end
                if pxx_max < max(th_pds)
                    axes_data(diseaseIndex,2) = max(th_pds);
                end
                if f_min > min(th_f)
                    axes_data(diseaseIndex,3) = min(th_f);
                end
                if f_max < max(th_f)
                    axes_data(diseaseIndex,4) = max(th_f);
                end
            end    
			
            %% store disease and frequency data
            eval([regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(dis_name,'[^\w'']',''),'_pds = th_pds;'])
            eval([regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(dis_name,'[^\w'']',''),'_f = th_f;'])
            
            %% save figure 
            if saveFig
                set(gcf,'PaperPositionMode','auto')
                print('-depsc2', '-loose',['.\Figures\multiDisease\',countryNames{countryIndex},'_',dis_name]);
                saveas(h, ['.\Figures\multiDisease\',countryNames{countryIndex},'_',dis_name,'.png']);
                close(h);
            end
						
            %% phase from fft
            Fs = 365/14; % sampling frequency, biweekly
            L  = length(caseDatai); % sample length
            N  = L ; %2^(nextpow2(L));
            Y  = fft(caseDatai, N);            
            NumUniquePts = ceil((N+1)/2);
            
            P2  = abs(Y/N);
            P1  = P2(1:NumUniquePts);

            if rem(N, 2) % odd nfft excludes Nyquist point
                P1(2:end) = P1(2:end)*2;
            else
                P1(2:end -1) = P1(2:end -1)*2;
            end            
            f = (0:(NumUniquePts-1))*Fs/N;
          
            figure(2)
            plot(f(2:end),P1(2:end),'bo-')
            xlabel('Frequency in years')
            ylabel(['FFT of ',dis_name,' in ', countryNames{countryIndex}]);
                
            x_candidates = find(abs(f-1)<1/12); % less than 1 month away from 1 year period
            if ~isempty(x_candidates)
                y_candidates = Y(x_candidates);
                phase = unwrap(angle(y_candidates));
                if length(y_candidates) > 1
                    fprintf('Multiple candidates are found. Frequency, PSD, and the corresponding phase are:\n')
                    for i=1:length(y_candidates)
                        fprintf('\t frequency %d: %g \t PSD: %g \t Phase: %g\n', i, f(x_candidates(i)), abs(y_candidates(i)), angle(y_candidates(i)));
                    end
                    fprintf('Phase of the max PSD is chosen\n')
                    y_max = find(y_candidates == max(abs(y_candidates)));
                    phase = phase(y_max);
                else
                    fprintf('\t frequency: %g \t PSD: %g \t Phase: %g\n', f(x_candidates), abs(y_candidates), angle(y_candidates));
                end
            else
                fprintf('No frequency near 1 year (one month apart) is found. Find the next best option\n');
                x_min = min(abs(f-1));
                x_candidates = find(abs(f-1)==x_min);
                y_candidates = Y(x_candidates);
                phase = unwrap(angle(y_candidates));
            end
            phaseInfo(countryIndex,diseaseIndex) = phase; % len_di/N is a correction factor for padding 0s
            if saveFig
                print('-depsc2', '-loose',['.\Figures\multiDisease\',countryNames{countryIndex},'_',dis_name,'_fft']);
                saveas(gcf, ['.\Figures\multiDisease\',countryNames{countryIndex},'_',dis_name,'_fft.fig']); 
            end
        end
        %% save disease and frequency data
        cd('IntermediateProcessingMatlabDataStructure')
        cd('CountryData')
        eval(['save ', regexprep(countryNames{countryIndex},'[^\w'']',''),'_disease_struct.mat ',regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind)]) 
        cd ..
        cd ..
    end
    cd('IntermediateProcessingMatlabDataStructure')
    save processed_info.mat chosenDN countryNames axes_data sth_ind 
    cd ..

end

%% processed data
cd('IntermediateProcessingMatlabDataStructure')
% load country names and disease names
load processed_info.mat
cd('CountryData')

for diseaseIndex=1:len_dis
    % create a struct for each disease
    dis_name = chosenDN{diseaseIndex};
    struct_name = [regexprep(dis_name,'[^\w'']',''),'_struct_',num2str(sth_ind)]; 
    eval([regexprep(dis_name,'[^\w'']',''),'_struct_',num2str(sth_ind),' = struct;'])
    
    for countryIndex=1:cnLen
        country_name = regexprep(countryNames{countryIndex},'[^\w'']','');
        eval(['load ',regexprep(countryNames{countryIndex},'[^\w'']',''),'_disease_struct'])
        eval([regexprep(dis_name,'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(countryNames{countryIndex},'[^\w'']',''),'_pds = ',regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(dis_name,'[^\w'']',''),'_pds;'])
        eval([regexprep(dis_name,'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(countryNames{countryIndex},'[^\w'']',''),'_f = ',regexprep(countryNames{countryIndex},'[^\w'']',''),'_struct_',num2str(sth_ind),'.',regexprep(dis_name,'[^\w'']',''),'_f;'])
    end
    eval(['save ', regexprep(dis_name,'[^\w'']',''),'_struct.mat ' regexprep(dis_name,'[^\w'']',''),'_struct_',num2str(sth_ind)])
end
cd ..
save phaseInfo 'phaseInfo'
save periodInfo 'periodInfo'
save sigInfo 'sigInfo'
close all



