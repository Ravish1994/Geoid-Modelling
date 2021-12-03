# Geoid Modelling Project CE678A

### 1. Study area:                                    --------------> [Done]
 Gravity for the Redefinition of the American Vertical Datum : Grav-D Project\
 Gravity data observed by the aircraft used during the project for Block_CN03\
 Block Extent: 40 to 43 degree North latitude && 96 to 102 degree West longitude\
 Region Name: Nebraska, USA\
 Data collection period: From 7-2014 to 9-2014\
 Data link: https://www.ngs.noaa.gov/GRAV-D/data_cn03.shtml

### 2. Reading the airborne gravity data of BLOCK CN03    --------------> [Done]
Raw_Data = importdata('CN03.txt')\
latitude_CN03 = Raw_Data.data(:,2)\
longitude_CN03 = Raw_Data.data(:,3)\
Ellht_cn03_NOAA = Raw_Data.data(:,4)\
gravity_CN03_NOAA = Raw_Data.data(:,5)

#### 2.1. Masking the airborne gravity data of BLOCK CN03 in the region of extent 40-41N and 97-98W
latnoaa = []\
lonnoaa = []\
g_NOAA = []\
h_aircraft1 = []\
for i = 1:length(latitude_CN03)\
    a = latitude_CN03(i)\
    b = longitude_CN03(i)\
    if (a>=40 & a<=41) & (b>=-98 & b<=-97)\
        latnoaa(end+1) = a\
        lonnoaa(end+1) = b\
        g_NOAA(end+1) = gravity_CN03_NOAA(i)\
        h_aircraft1(end+1) = Ellht_cn03_NOAA(i)\
    end\
end

#### 2.2. Creating equispaced latitude longitude points using linspace function of NOAA data of CN03 block
[X_lat Y_lon] = meshgrid(double(linspace(40,41,100)),double(linspace(-98,-97,100)))

### 3. Calculating DEM for NOAA data observation points       --------------> [Done]
SRTM DEM data:\
Resolution: 1 arc-second for global coverage (~30 meters)\
Extent: 40 to 43 degree North latitude && 96 to 102 degree West longitude\
Region name: Nebraska, USA\
Data collection period: 2014\
Number of tiles available: 40\
Data link: https://earthexplorer.usgs.gov/

### 4. Combining all the raw DEM tiles for getting DEM over throughout the extent of block CN03   --------------> [Done]
DEM_raw_CN03 = [v_8 v_7 v_6 v_5 v_4 v_3 v_2 v_1;...\
                u1_8 u1_7 u1_6 u1_5 u1_4 u1_3 u1_2 u1_1;...\
                v8 v7 v6 v5 v4 v3 v2 v1;...\
                u8 u7 u6 u5 u4 u3 u2 u1;...\
                t8 t7 t6 t5 t4 t3 t2 t1]\

#### 4.1. Reading Geotiff file using geotiffread function\
F11 = 'n40_w098_1arc_v3.tif';\
t11 = geotiffread(F11);  It provides integer value of the DEM data \
Here we tried by readGeoraster function too It gives same integer values\
Since while interpolating the DEM data we need double values\ 
T11 = double(t11);\
Since DEM data is of size 3601x3601 so we need to create the grid of\
co-ordinates corresponding to those DEM points\

#### 4.2. Creating meshgrid of SRTM DEM coordinates
a1 = double(linspace(40,41,3601));    % Since DEM Extent is 40-41 North Latitude\
a2 = double(linspace(-98,-97,3601));  % Since DEM Extent is 97-98 West Longitude\
[Lat_DEM,Lon_DEM] = ndgrid(a1,a2);    % Creating Meshgrid of Coordinates\

#### 4.3. Interpolating DEM data over the latitude and longitude of the NOAA 
Since Observation points masked for the DEM extent comes inside the grid of DEM coordinates so griddata function will give best\ interpolation results on those observation points , For the best result we will use cubic method of interpolation as it is the best one in\ all other methods like 'nearest', 'natural' and 'v4'.\

h_DEM_NOAA = double(griddata(Lon_DEM,Lat_DEM,T11,lonnoaa,latnoaa,'cubic'));\
Interpolating over created meshgrid of Latitude and Longitude of airborne gravity data\
h_DEM_NOAA_mesh = double(griddata(Lat_DEM,Lon_DEM,T11,X_lat,Y_lon,'cubic'));\

### 5. Calculation of Free Air Correstion(FAC) Over the latitude and the longitude of gravity Observation points of the NOAA  --------------> [Done]
Interpolating the ellipsoidal height of the observed gravity data over created latitude and longitude meshgrids of NOAA data sets\
H_aircraft = double(griddata(latnoaa,lonnoaa,h_aircraft1,X_lat,Y_lon,'v4'));\

Finding the height of aircraft over topographical surface Since we need to provide correction for the air between the topography and the\ measurement points We can find this substracting the DEM from the ellipsoidal height\
H_air = H_aircraft-h_DEM_NOAA_mesh;\

Free air correction can be computed as:\
FAC_NOAA = 0.3086*H_air;

### 6. Calculation of Geoid Undulations Over the block CN03 from Geoid undulations of EGM2008  --------------> [Done]
Reading the EGM008 file\
b = 'EGM2008_N.txt';\
A1 = importdata(b);\
Lon_EGM08 = A1(:,1)-360; % Since Longitude is more than 180 degree\
Lat_EGM08 = A1(:,2);\
N_EGM08  = A1(:,3);\

Interpolating the Geoid undulation value of EGM2008 over meshgrid of Latitude and Longitude of airborne gravity data\
N_NOAA_mesh = double(griddata(Lon_EGM08,Lat_EGM08,N_EGM08,Y_lon,X_lat,'cubic'));\

### 7. Calculation of Bouguer Anomaly over the block CN03   --------------> [Done]
Height for calculating Bouguer anomaly\
H_bg = h_DEM_NOAA_mesh - N_NOAA_mesh;\
Calculating Bouguer Anomaly in mGal using standard formula\
G_NOAA_mesh = double(griddata(lonnoaa,latnoaa,g_NOAA,Y_lon,X_lat,'v4'));\
Bg_NOAA_mesh = G_NOAA_mesh - 0.1119*H_bg + FAC_NOAA;\
 
### 8. Calculating normal gravity value(gamma) in mGal on CN03 block  --------------> [Done]
WGS84 Parameters for calculating the normal gravity value(gamma)\
a = 6378137 ;             % In meter\ 
b = 6356752.3142 ;        % In meter\
g_equator = 983218.49378; % Normal gravity value at the equatorial point\
g_pole = 978032.7;        % Normal gravity value at the polar point\

### 9. Calculating the normal gravity value(gamma)     --------------> [Done]
g11 = (a*g_equator*(cosd(X_lat)).^2+b*g_pole*(sind(X_lat)).^2);\
g21 = sqrt(a^2*(cosd(X_lat)).^2+b^2*(sind(X_lat)).^2);\
gama = g11./g21;\
R = radius of earth;\
C = R./(2*gama1);\

### 10. Calculating the gravity disturbances  --------------> [Done]
T_NOAA_mesh = G_NOAA_mesh - gama;\

### 11. Calculating the free air gravity anomalies in mGal on block CN03  --------------> [Done]
FDg_NOAA_mesh = G_NOAA_mesh + FAC_NOAA - gama;\

### 12. Adding bundles for the spherical harmonic coefficients Cnm and Ylm calculations --------------> [Done]
addpath(genpath('D:\Physical_Geodesy\Geoid_Modelling\Software_Ellman_2004\SHbundle-master'));\
addpath(genpath('D:\Physical_Geodesy\Labs\Lab_05\uberall-master'));\
addpath(genpath('D:\Physical_Geodesy\Labs\Lab_05\gracebundle-master'));\

### 13. Finding the gravity anomaly values delG_N(EGM)  --------------> [Done]
Importing GGM data for degree and order and Cnm\
load GGM01S.prn\
mat = GGM01S(:,1:4);\

Using standard formula calculated as in the lab document we get:\
As we increase the degree gravity anomaly will decrease as we can see by changing the value of n Therefore we need to take the\
maximum degree to get the less gravity anomaly values:  n = 120 Gravity Anomaly is the value which is difference between our observed gravity values and predicted by the model GGM01S\

delgn_NOAA = [];\
for i = 1:length(lonnoaa)                           % Data range\
    f_nm = 0;\
    n = 120;                                        % Maximum degree given in GGM data\
    for m = 0:n                                   % Since order varies from 0 to n\
        a  = 6378137;                               % Semimajor axis of the WGS84 reference ellipsoid\
        gm = 3986005*10^8;                          % Geocentric gravitational constant of the earth\
        r  = 6400000;                               % Geocentric radius of the earth\
        C1 = (gm/(a.^2));\
        C2 = (a/r).^(n+2);\
        C3 = (n-1);\
        Ct2  = C1*C2*C3;\

        % Cnm : Fully normalised spherical harmonic coefficient of\
        % disturbing potential from the normal gravity field refering to\
        % the radius of a\

        C_nm = get_coeff(mat,n,m);                  % Getting the Cnm coefficient from GGM files\
        P_nm = plm(n,m,latnoaa(i)*pi/180);          % Calculating Associated legendre polynomials\

        % Ynm : Fully normalised surface spherical harmonic\

        Y_nm  = P_nm*cosd(m*lonnoaa(i)*pi/180);\   
        f_nm  = f_nm + Ct2*C_nm*Y_nm;\
    end\
    delgn_NOAA(end+1)= f_nm;\
end\
dgN = double(griddata(latnoaa,lonnoaa,delgn_NOAA,X_lat,Y_lon,'cubic'));\

### 14. Finding Q_NLsi : Molodenskii’s truncation coefficient   --------------> [Done]
Input data \
PSI, Integration radius for geoid determination, 0.2\
Max_exp,Maximum degree of series expansion, e.g. M_max=121\
L_max, Maximum modification degree, e.g. L_max = 70\

Output data\
ENK, Paul's coeff. are presented in the matrix with dimens. (Max_exp-1)x(L_max-1) \
QN, Molodenskii's trunc. coeff. are presented as a vector of (Max_exp-1) elements\
Passing 120 as maximum degree \
[ENK,QN]=trunc_coeff_3(0.2,120,70);\

### 15. Finding modified Stoke's Kernel(SLG_bias)    --------------> [Done]
Input data\

PSI, Integration radius for geoid determination\
L_max, Maximum modification degree, e.g. L_max = 70\
M_max, Upper limit of the used geopotential coefficients, e.g. M_max = 70\
Sn_bias, Vector of biased LS modification parameters\
Sn_k_unb, Matrix of unbiased LS modification parameters (columnwise, solved from T-SVD)\
Sn_k1_unb, Matrix of unbiased LS modification parameters (columnwise, solved from T-TLS)\
Sn_k_opt, Matrix of optimum LS modification parameters (columnwise, solved from T-SVD)\
Sn_k1_opt, Matrix of optimum LS modification parameters (columnwise, solved from T-TLS)\

For Sn and Bn coefficients use software by Ellmann ''LS_coeff_v2.m''\
Sn_bias = importdata('Sn_bias.prn');\
Sn_k_unb = importdata('Sn_k_unb.prn');\
Sn_k_opt = importdata('Sn_k_opt.prn');\
Sn_k1_unb = importdata('Sn_k1_unb.prn');\
Sn_k1_opt = importdata('Sn_k1_opt.prn');\
Bn_opt = importdata('Bn_opt.prn');\
PSI = 0.2;\
M_max =120;\
L_max = M_max;\
[SL,SLG_bias] = modif_STOKES_kernel(PSI,M_max,L_max,Sn_bias,Sn_k_unb, Sn_k_opt,Sn_k1_unb,Sn_k1_opt);\

### 16. Finding the approximate Geoid Undulation value N   --------------> [Done]
N_Bare = [];\
for i = 1:length(lonnoaa)\
    sum = 0;\
    for n = 2:120\
        N1 = C(i)*((2/(n-1))-QN(n-1)-Sn_k_opt(n-1)+Bn_opt(n-1))*dg_N(i);\
        sum = sum+N1; %% Converting g_N in mGal\
    end\
    N_Bare(end+1) = sum;\
end\
N_Bar = griddata(latnoaa,lonnoaa,N_Bare,X_lat,Y_lon,'nearest');\

### 17. Providing Additive corrections   

#### 17.1. Topographic Combined correction  --------------> [Done]
Rho = 2.65; % Crust density of earth in g/cc\
It is recommended to take Rho in g/cc\
gama in mGal\
H_bg : Orthometric height\
dN_topo_Combined = ((-2*pi*(6.67*10^(-11))*Rho*H_bg^2)./gama)*10^5;\

#### 17.2. Ellipsoidal correction   --------------> [Done]
Integration cap PSI in radian\
si_0 = 0.2;\
N_bar = Approximate geoid height\
dN_Ellipsoidal = (si_0)*((0.12-0.38*(cosd(90-X_lat)).^2).*FDg_NOAA_mesh+0.17*N_Bar.*(sind(90-X_lat))^2);\

#### 17.3. Atmospheric Combined Correction --------------> [Done]
dN_atm_Combined = [];\
for i = 1:length(lonnoaa)\
    sum = 0;\
    for n = 2:120\
        rho_0 = 1.23;\
        Cg = (2*pi*r*rho_0/(gama1(i)));\
        C2 = QN(n-1);\
        CS = Sn_k1_opt(n-1);\
        % As Refered to Abdalla 2011 paper to minimise the error they took\
        % H_np = 0.00002 meter\
        H_nP = 0.00002;\
        C3 = 3*(n+1)/(2*n+1);\
        C4 = 4/(n-1);\
        Ceff = C3*C2+CS-C4;\
        DN = Cg*Ceff*H_nP;\
        sum = sum+DN;\
    end\
    dN_atm_Combined(end+1)= sum/10^5;\
end\
N_atm_Combined = double(griddata(latnoaa,lonnoaa,dN_atm_Combined,X_lat,Y_lon,'nearest'));\

#### 17.4. Downward Continuation Correction --------------> [Done]
DNW_Correct = [];\
for i = 1:length(lonnoaa)\
    sum1 = 0;\
    for n = 2:120\
        C1 = r/(2*gama1(i));\
        C2 = QN(n-1);\
        C3 = Sn_k1_opt(n-1);\
        C4 = (r/(r+h_DEM_NOAA(i)))^(n+2)-1;\
        C5 = delgn_NOAA(i);\
        DN = C1*(C2+C3)*C4*C5;\
        sum1 = sum1+DN;\
    end\
    DNW_Correct(end+1)= sum1;\
end\
DNW_Correction = double(griddata(latnoaa,lonnoaa,DNW_Correct,X_lat,Y_lon,'nearest'));\

#### 17.5. Estimated Geoid Model   --------------> [Done]
N_cap = N_Bar + dN_topo_Combined + dN_Ellipsoidal + N_atm_Combined +DNW_Correction;\
xlabel('Longitude')\
ylabel('latitude')
caxis([min(min(N_NOAA_mesh)),max(max(N_NOAA_mesh))])
| STEPS | Progress |
| ------ | ------ |
| Study Area| Done|
| Data Visulisation | Done |
| Computation of LSMSK parameters from Ellamann software| Done|
| DEM Height computation| Done|
| Gravity Anomaly(Free Air Bouguer)| Done |
| Geoid Undulation(Approx)| Done |
|  Additive Correction| Done |
| Geoid estimated| Done |


