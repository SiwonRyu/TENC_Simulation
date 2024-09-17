clear all
baseroot = 'C:\rsw\Replication_CP_2021\Empirical Illustration\';
Dyadic_data_raw = importdata([baseroot,'Dyadic_data.xlsx']);
Ind_data_raw = importdata([baseroot,'Ind_data.xlsx']);

Data_dyadic = Dyadic_data_raw.data;
Data_individual = Ind_data_raw.data;
dyadic_index = Data_dyadic(:,2)~=Data_dyadic(:,3);
Data_dyadic = Data_dyadic(dyadic_index,:);



%%

[est_zeta, est_beta, est_beta_L, est_pi] = est_RE(10,Data_dyadic,Data_individual)

