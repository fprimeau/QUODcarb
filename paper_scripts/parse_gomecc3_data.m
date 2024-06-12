
% Download GOMECC-3 Data from website and parse into data structure

% download link found at this website:
% https://www.aoml.noaa.gov/ocd/gcc/GOMECC3/
url = 'https://www.aoml.noaa.gov/ocd/gcc/GOMECC3/GOMECC-3%20Cruise%20Data.xlsx';

T = webread(url); % reads the table into Matlab

nD = length(T.Sample_ID);

% create empty vectors to receive the data
depth = [];
pres = [];
temp = [];
sal = [];
TC = [];
TA = [];
co3 = [];
ph_tot = [];
pco2 = [];
TP = [];
TSi = [];

% check if each datapoint has all the measurements and parse into vectors
% flags: 2 = good, 6 = mean of replicates
for i = 1:nD
    if T.CARBONATE_FLAG_W(i) == 2 || T.CARBONATE_FLAG_W(i) == 6
        if T.PH_TOT_FLAG_W(i) == 2 || T.PH_TOT_FLAG_W(i) == 6
            if T.PCO2_FLAG_W(i) == 2 || T.PCO2_FLAG_W(i) == 6
                if T.PHOSPHATE_FLAG_W(i) == 2 && T.PHOSPHATE_UMOL_KG(i) >= 0 ...
                        || T.PHOSPHATE_FLAG_W(i) == 6 && T.PHOSPHATE_UMOL_KG(i) >= 0
                    if T.SILICATE_FLAG_W(i) == 2 || T.SILICATE_FLAG_W(i) == 6
                        if T.TA_FLAG_W(i) == 2 || T.TA_FLAG_W(i) == 6
                            if T.DIC_FLAG_W(i) == 2 || T.DIC_FLAG_W(i) == 6
                                co3(end+1) = T.CARBONATE_UMOL_KG(i);
                                ph_tot(end+1) = T.PH_TOT_MEA(i);
                                pco2(end+1) = T.PCO2_MEA_UATM(i);
                                TP(end+1) = T.PHOSPHATE_UMOL_KG(i);
                                TSi(end+1) = T.SILICATE_UMOL_KG(i);

                                depth(end+1) = T.DEPTH_METER(i);
                                pres(end+1) = T.CTDPRS_DBAR(i);
                                temp(end+1) = T.CTDTMP_ITS90_DEG_C(i);
                                sal(end+1) = T.CTDSAL_PSS78(i);
                                TC(end+1) = T.DIC_UMOL_KG(i);
                                TA(end+1) = T.TA_UMOL_KG(i);
                            end
                        end
                    end
                end
            end
        end
    end
end

% parse into a single matrix
data(1,:) = sal;
data(2,:) = temp;
data(3,:) = pres;
data(4,:) = depth;
data(5,:) = TC;
data(6,:) = TA;
data(7,:) = TP;
data(8,:) = TSi;
data(9,:) = ph_tot;
data(10,:) = pco2;
data(11,:) = co3;

% save for future use:
save data.mat data;

% note three of the phosphate measurements were negative, cut them out
