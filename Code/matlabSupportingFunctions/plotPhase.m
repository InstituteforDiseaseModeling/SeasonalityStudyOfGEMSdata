% plot phase info
% Min Roh

close all
clear
clc
cd('DataStructure')
load('phaseInfo.mat')
load('sigInfo.mat')
load('periodInfo.mat')
cd ..

chosenInds = [0, 2, 3, 5, 6, 8, 9, 10]+4;
chosenDN   = {'Eligible', 'Campylobacter', 'Shigella', 'ETEC ST', 'tEPEC', 'Cryptosporidium',  'Rotavirus', 'Adenovirus'};
len_dis    = length(chosenInds);


%% determine exact peak and label
% compute the start date of each disease
data              = readtable('../Data/estimated-positive-fortnightly.csv');
countryNames      = unique(data.country);
% reorder countries west to east: The Gambia, Mali, Kenya, Mozambique, Pakistan, India, Bangaldesh)
countryInds       = [7 4 3 5 6 2 1];
[m,n]             = size(phaseInfo);
start_dates       = zeros(m,n);
for i=1:m
    d_start          = find(strcmp(data.country,countryNames{i}),1);
    start_dates(i,:) = datenum(data{d_start,3});
end

%% phase translation, value set to radian of the peak
m_format = 'mm';
d_format = 'dd';
for i=1:m
    for j=1:n
        m_ij = str2double(datestr(start_dates(i,j),m_format));
        d_ij = str2double(datestr(start_dates(i,j),d_format));       
        md_radian = 2*pi - ((m_ij-1)*30.417+d_ij) * 2*pi / 365;
        
        p_ij = phaseInfo(i,j);
        if p_ij <0            
            phaseInfo(i,j) = rem(phaseInfo(i,j)*-1 - md_radian+2*pi, 2*pi);
        else
            phaseInfo(i,j) = rem((phaseInfo(i,j)-2*pi)*-1 - md_radian+2*pi, 2*pi);
        end
    end
end

%% rearrange
phaseInfo   = phaseInfo(countryInds,:)';
start_dates = start_dates(countryInds,:)';
sigInfo     = sigInfo(countryInds,:)';
periodInfo  = periodInfo(countryInds,:)';

%% superimpose a dot whose relative size indicates strong annual seasonality
col_mat = [0 1 1;
          .31 .31 .31];
m_size = 20;      
f_size = 12;

%% convert to numbers, whole day
phaseInfo_datenum = round(phaseInfo*365/(2*pi));
formatOut         = 'mm-dd';
phaseInfo_datestr = datestr(phaseInfo_datenum,formatOut);


%% plot phase
figure(3)
imagesc(phaseInfo)
alpha(0.8)
xticklabels(countryNames(countryInds));
yticks(1:length(chosenDN))
yticklabels(chosenDN);
set(gca,'TickLabelInterpreter','none');
set(gca,'fontweight','bold','fontsize',12);
set(gcf, 'Position', [50 120 1100 630])
hold on
% custom colormap for circular map
slen = 200;
cmat = [1 .333 0; % red
    .988 .812 .188;  % yellow 
    .29 .796 .518; % green    
    0.2745 0.231 0.871; % blue
    .667 .333 1]; % purple

seg1 = [linspace(cmat(1,1),cmat(2,1),slen) linspace(cmat(2,1),cmat(3,1),slen) linspace(cmat(3,1),cmat(4,1),slen) linspace(cmat(4,1),cmat(5,1),slen) linspace(cmat(5,1),cmat(1,1),slen)]'; %linspace(cmat(6,1),cmat(1,1),slen)]';
seg2 = [linspace(cmat(1,2),cmat(2,2),slen) linspace(cmat(2,2),cmat(3,2),slen) linspace(cmat(3,2),cmat(4,2),slen) linspace(cmat(4,2),cmat(5,2),slen) linspace(cmat(5,2),cmat(1,2),slen)]'; %linspace(cmat(6,2),cmat(1,2),slen)]';
seg3 = [linspace(cmat(1,3),cmat(2,3),slen) linspace(cmat(2,3),cmat(3,3),slen) linspace(cmat(3,3),cmat(4,3),slen) linspace(cmat(4,3),cmat(5,3),slen) linspace(cmat(5,3),cmat(1,3),slen)]'; %linspace(cmat(6,3),cmat(1,3),slen)]';
cmap = [seg1 seg2 seg3]; 
colormap(cmap)
caxis([0 2*pi])
months = linspace(0, 2*pi, 13);
cb = colorbar('Ticks',months,...
    'TickLabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan'});

%% write the dates
% 1 - not significant
% 2 - 10%
% 3 - 5%
% 4 - 1%
% 5 - 0.1%
for i=1:m
    for j=1:n
        m_ij = sigInfo(j,i);
        
        %         if m_ij == 2
        %             text(i-.06,j-.2,'*','Color',col_mat(1,:),'FontSize',m_size)
        %         elseif m_ij == 3
        %             text(i-.06,j-.2,'**','Color',col_mat(1,:),'FontSize',m_size)
        %         elseif m_ij == 4
        %             text(i-.11,j-.2,'***','Color',col_mat(1,:),'FontSize',m_size)
        %         elseif m_ij == 5
        %             text(i-.15,j-.2,'****','Color',col_mat(1,:),'FontSize',m_size)
        %         end
%         if m_ij == 2 || m_ij == 3
%             text(i-.06,j-.2,'*','Color',col_mat(1,:),'FontSize',m_size)
%         elseif m_ij == 4
%             text(i-.11,j-.2,'**','Color',col_mat(1,:),'FontSize',m_size)
%         elseif m_ij == 5            
%             text(i-.15,j-.2,'***','Color',col_mat(1,:),'FontSize',m_size)
%         end
        text(i-.18,j+.05, datestr(phaseInfo_datenum(j,i),formatOut),'FontSize',f_size) %, '(',num2str(periodInfo(j,i)*365-365,'%2.f'),')'],'FontSize',f_size) 
    end
end
% plot(-1:8,linspace(8.5,8.5,10),'--b','LineWidth',2)
% plot(-1:8,linspace(6.5,6.5,10),'--g','LineWidth',2)
% plot(-1:8,linspace(1.5,1.5,10),'--r','LineWidth',2)
% save figure
set(gcf,'PaperOrientation','landscape');
print('-depsc2', '-loose', '.\Figures\phase\phaseInfo_west2east_4sigs');
saveas(gcf, '.\Figures\phase\phaseInfo_west2east_4sigs.pdf');
hold off