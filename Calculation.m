%% Loading Data

load SCEDC_1982_to_2022.mat

%% Declustering 

SCEDC_raw_d = Gardner_Knopoff(SCEDC_raw,mainshock,3);
%^removes foreshocks and aftershocks of M7.0+ earthquakes 

%% Monthly b_ML, b_LS, D2 Calculations (with declustering) 

bot = 2.5;
top = 4.5;

SCEDC_temp_d = SCEDC_raw_d(SCEDC_raw_d(:,7) >= bot & SCEDC_raw_d(:,7) <= top,:);
%^restricts data to earthquakes with magnitudes contained in the interval
%[bot, top]

month_count = (max(SCEDC_temp_d(:,3)) - min(SCEDC_temp_d(:,3)) + 1)*12;

b_ML_d = NaN(1,month_count);
b_ML_error_d = NaN(1,month_count);

b_LS_d = NaN(1,month_count);
b_LS_error_d = NaN(1,month_count);

D2_d = NaN(1,month_count);
D2_error_d = NaN(1,month_count);

s = 10.^(-1:0.1:3);
i = 11:1:21;
%^s(i) corresponds to s-values contained in the interval [1 km, 10 km]

count = 1;

for ii = min(SCEDC_temp_d(:,3)):1:max(SCEDC_temp_d(:,3))
    year = SCEDC_temp_d(SCEDC_temp_d(:,3) == ii,:);
    for jj = 1:1:12
        month = year(year(:,1) == jj,:);
        if size(month,1) >= 30
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % b_ML value calculation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            b_ML_d(count) = b_calc_fun(month(:,7));
            b_ML_error_d(count) = b_error_calc_fun(month(:,7),b_ML_d(count));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % b_LS value calculation 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            magnitudes = month(:,7);
            
            Fanta = NaN(length(unique(magnitudes)),2);
            Fanta(:,2) = unique(magnitudes);
            %^magnitudes are arranged from top to bottom in ascending order

            for mm = 1:1:length(unique(magnitudes))
                Fanta(mm,1) = sum(magnitudes >= Fanta(mm,2));
            end

            b_LS_linear = polyfit(Fanta(:,2),log10(Fanta(:,1)),1);
            b_LS_d(count) = -b_LS_linear(1);

            Sunkist = Fanta';

            b_LS_error_d(count) = error_calc_fun(Sunkist(2,:),log10(Fanta(:,1)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % D2 value calculation 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            CI = CI_calc_fun(month(:,8:10),s);

            D2_d(count) = D2_calc_fun(CI,s,i);
            D2_error_d(count) = error_calc_fun(log10(s(i)),log10(CI(i)));

        end
        if count ~= month_count
            count = count + 1;
        end
    end
end

%% Saving Data 

filename = 'Temporal_Results_After_Declustering.mat';
save(filename)