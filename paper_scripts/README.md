## Fennell & Primeau 2024
To replicate the figures in Fennell & Primeau 2024, first run `parse_gomecc3_data.m`. This will download the raw datafile to your local directory and parse it into `data.mat`. Next, run `driver.m` to run through the GOMECC-3 dataset (Barbero et al., 2019) parsed into `data.mat` 26 times, one for each possible combination. Please create a folder named `output_mat_files` to contain the outputs. The file `driver.m` takes 20GB and 4.5 hours on a remote computer for M. Fennell. 

If you are only interested in figure 4, figure 5, and/or figure B7 then you can stop the driver after line 81 (starts with `save output_mat_files/est26.mat`). Figures 2, 3, and 6 require more than just the fully over-determined run which is saved as `est26.mat`. 


Barbero, L., Pierrot, D., Wanninkhof, R., Baringer, M.O., Hooper, J.A., Zhang, J.Z., Smith, R.H., Byrne, R.J., Langdon, C., Herna ́ndez-Ayo ́n, M., Schnetzer, A., Stauffer, B., Herzka, S., Cano-Compair ́e, J., Hern ́andez, F.J., Pech, D., Hu, C., English, D., Ondrusek, M., Olascoaga, J., Derr, G., 2019. Third Gulf of Mexico Ecosystems and Carbon Cycle (GOMECC-3) Cruise URL: https://repository.library.noaa.gov/view/noaa/ 21256, doi:10.25923/Y6M9-FY08. publisher: Atlantic Oceanographic and Meteorological Laboratory.



