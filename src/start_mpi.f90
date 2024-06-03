!*************************************************************************************************
! Start up procedure for multiple grid cells
! using parallel computation with MPI
!*************************************************************************************************
PROGRAM omp_kickoff
   USE mpi
   USE data_structure
   implicit none
   
!_____________ Set Parameters
!Directory of spin up files
   character(len=*),parameter::Loc_spnin_files='/data/hsato/siberia/spnin1850/'

!Directory of aCO2 data
   character(len=*),parameter::Loc_CO2_files='co2_holocene.dat' !RCP8.5
!  character(len=*),parameter::Loc_CO2_files='co2_1850_2100_rcp26.dat' !RCP2.6

!Directory of climate data for future projection
 character(len=*),parameter::Loc_climate_data2='/homedata/ibertrix/Climate_data/EWEMBI-Europe-SEIB_1979-2016/'

 !character(len=*),parameter::Loc_climate_data2='/work/G10203/hsato/climate_MirocAR5Base_daily_0.5deg_V3/RCP85/'
 !character(len=*),parameter::Loc_climate_data2='/work/G10203/hsato/climate_MirocAR5Base_daily_0.5deg_V3/RCP26/'

!Directory of climate data for histrical reconstruction
! character(len=*),parameter::Loc_climate_data1='/work/G10203/hsato/climate_MirocAR5Base_daily_0.5deg_V3/historical/'
 character(len=*),parameter::Loc_climate_data1='/homedata/ibertrix/Climate_data/EWEMBI-Europe-SEIB_1979-2016/'   

!Set year number of climate & CO2 data
   integer,parameter::YearMaxClimate1 = 38  !historical(1850-2005)
   integer,parameter::YearMaxClimate2 =  0  !RCPs (2006-2100)
   integer,parameter::YearMaxCO2      = 300  !Year length of the CO2 data (1850-2100)
   
!Land mask file
   character(len=*),parameter::Fn_landmask = &
                               'landmask_0.25deg.txt'
   
!Directory for writing an output file of each simulation grid
   character(len=*),parameter::Loc_result_files = &
                               './result_tmp/'
   
!Directory for writing analysis data for whole simulation grid
   character(len=*),parameter::Loc_analysis_files= &
                               './result_visualize/'
   
!Location of soil property data
!A same data-set is employed at every simulation grids
!This data-set was obtained from GSWP2 (http://www.iges.org/gswp2/)
   character(len=*),parameter::Fn_landprop = 'land_prop.txt'
   
!Maximum grid number for longitude and latitude
  !integer,parameter::LatMax = 360 !Grid number for latitude  (@ 0.5deg Grid System)
  !integer,parameter::LonMax = 720 !Grid number for longitude (@ 0.5deg Grid System)
  integer,parameter::LatMax = 720 !Grid number for latitude  (@ 0.25deg Grid System)
  integer,parameter::LonMax = 1440 !Grid number for longitude (@ 0.25deg Grid System)
   
!Grid length for simulation
  !(East Siberia @ 0.5deg grid mesh)
  integer,parameter::LatNoStart =  73 !(N55)
  integer,parameter::LatNoEnd   = 228 !(N60)
  integer,parameter::LonNoStart = 665 !(E120)
  integer,parameter::LonNoEnd   = 920 !(E135)
  
!_____________ Set Variables
   !MPI control variables
   integer myid, numprocs, ierr
   
   !Coodinate variables for parallel computation
   integer gridNo, gridNoMax        !Sequential number for simulation grid cell
   integer latNo, lonNo             !Sequential number for north-south and west-east
   real    lat,lon                  !Latitude and Longitude
   integer,dimension(LatMax*LonMax)::latNo_save, lonNo_save
   
   !For preparing Land-Ocean mask
   logical,dimension(LonMax,LatMax)::Landmask     !landmask; true->land, false->water
   integer,dimension(1:LonMax)     ::readerLonMax !For reading land mask
   integer,dimension(360*180)::Mask, SoilClass
   real   ,dimension(360*180)::ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt
   integer point
   
   !Atomospheric CO2 time-series @ ppm
   real,dimension(300)::aco2_1850to2100
   integer,parameter::startyear = 300
   integer,parameter::co2start = 12000 - startyear
   integer,parameter::co2stop = co2start + YearMaxCO2 - 1
   !Counter
   integer i
   
!_____________ Read Parameters
!Read Parameter files
   open (1, file='parameter.txt', action='READ', status='OLD')
      read ( unit=1, nml=Control)       
      read ( unit=1, nml=PFT_type)      
      read ( unit=1, nml=Respiration)   
      read ( unit=1, nml=Turnover_n)    
      read ( unit=1, nml=Metabolic)     
      read ( unit=1, nml=Assimilation)  
      read ( unit=1, nml=Dynamics)      
      read ( unit=1, nml=Disturbance)   
      read ( unit=1, nml=Soil_resp)     
   close (1)
   
!______________ Prepare landmask (extension for wide area simulation)
!Read Landmask
   landmask(:,:) = .false.
   !landmask(:,:) = .true.
   open (1, file=Fn_landmask, status='OLD')
   do latNo=1, LatMax
      read(1,*) readerLonMax(1:LonMax)
      do lonNo=1, LonMax
         if (readerLonMax(lonNo)==1) landmask(lonNo,latNo)=.true.
      end do
   end do
   close(1)
   
   Mask        (:)=0
   ALT         (:)=0.0
   Albedo_soil0(:)=0.0
   W_sat       (:)=0.0
   W_fi        (:)=0.0
   W_mat       (:)=0.0
   W_wilt      (:)=0.0
   SoilClass   (:)=0
   open (1, file='land_prop.txt', status='OLD')
   do i=1, 180*360
      read(1,*) Mask(i), ALT(i), Albedo_soil0(i), &
                W_sat(i), W_fi(i), W_mat(i), W_wilt(i), SoilClass(i)
   end do
   close(1)
   
   !Landmask correction
   do latNo=1, LatMax
   do lonNo=1, LonMax
      if (.not. landmask(lonNo,latNo)) cycle
      
      lat   =    90.0 - (real(latNo)-0.5) * 0.25
      lon   = - 180.0 + (real(lonNo)-0.5) * 0.25 !180.0 + (real(lonNo) - 0.5) * grid_spacing
      point = (90-int(lat)-1)*360 + int(lon+180) + 1 !grid point number @1.0 deg system
      
      if (Mask     (point)==0) landmask(lonNo,latNo)=.false.
      if (SoilClass(point)==0) landmask(lonNo,latNo)=.false.
      if (SoilClass(point)==9) landmask(lonNo,latNo)=.false.
   end do
   end do
   
!______________ Read time-series of atmospheric CO2 concentration
   Open (1, file=trim(Loc_CO2_files), status='OLD')
   do i = 1, co2start-1
        read(1, *) aco2_1850to2100(1)
   end do
   do i = co2start, co2stop !1850~2100
      read(1, *) aco2_1850to2100(i-co2start+1)
   end do
   Close (1)
   
!______________ Provide reference IDs for each grid cell
   gridNo = 0
   do latNo = LatNoStart, LatNoEnd
   do lonNo = LonNoStart, LonNoEnd
      if (.not. landmask(lonNo,latNo)) cycle
      
      gridNo             = gridNo + 1
      gridNoMax          = gridNo
      latNo_save(gridNo) = latNo
      lonNo_save(gridNo) = lonNo
   end do
   end do
   if (gridNoMax==0) Stop
   
!______________ Initialization of MPI
   Call MPI_INIT( ierr )
   Call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
   Call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
   
!______________ Conduct parallel computation with MPI
   if (myid==1) write(*,*) 'Start parallel computation....'
   
   DO gridNo = 1, gridNoMax
   if ( myid == mod(gridNo-1, numprocs) ) then
      Call start(myid, latNo_save(gridNo), lonNo_save(gridNo), aco2_1850to2100, &
                 YearMaxClimate1, YearMaxClimate2, YearMaxCO2, &
                 Loc_climate_data1, Loc_climate_data2, Loc_result_files, Loc_spnin_files)
   endif
   END DO
   
!______________ Termination procudure of MPI compuation
   Write(*,*) 'Terminate processor:', myid, ', Total processor:', numprocs
   Call MPI_Barrier(MPI_COMM_WORLD, ierr) !Wait until end of all simulation
   Call MPI_FINALIZE(ierr)                !Finalize procudure for MPI
   
!_____________ Conversion of output files for each simulation grids into maps
!ObhÊÌt@CoÍðA±±Ån}ÉÁH·é
   if (myid==0) then
      write(*,*) 'Converting result files....'
      !Analysis output values of Last 10 yers
      Call after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   endif
   
STOP
END PROGRAM omp_kickoff



!*************************************************************************************************
! Start up procedure for each simulation grid
!*************************************************************************************************
Subroutine start (myid, latNo, lonNo, aco2_1850to2100, &
                  YearMaxClimate1, YearMaxClimate2, YearMaxCO2, &
                  Loc_climate_data1, Loc_climate_data2, &
                  Loc_result_files, Loc_spnin_files)
   
   USE data_structure
   implicit none
   
!_____________ Set Augment
   integer            ,intent(IN):: myid, latNo, lonNo
   integer            ,intent(IN):: YearMaxClimate1, YearMaxClimate2, YearMaxCO2
   character(len=*)   ,intent(IN):: Loc_climate_data1, Loc_climate_data2
   character(len=*)   ,intent(IN):: Loc_result_files, Loc_spnin_files
   real,dimension(YearMaxCO2),intent(IN):: aco2_1850to2100 !Atomospheric co2 concentration (ppm)
   
!_____________ Set Variables
!Climate data
   real,allocatable,dimension(:,:)  ::&
    tmp_air       ,& !1. Surface air temperature (Celcius)
    tmp_air_range ,& !2. Daily range of tmp_air (Celcius)
    prec          ,& !3. Precipitation (mm day-1)
    rad_short     ,& !4. Shortwave radiation, downward @ midday (W m-2)
    rad_long      ,& !5. Daily mean of longwave radiation, downward (W m-2)
    wind          ,& !6. Wind velocity (m s-1)
    rh               !7. Relative humidity (%)
   
   real,allocatable,dimension(:,:,:)::&
    tmp_soil      !Soil temperature for each layers (Celcius)
   
!Location data
   integer::Mask         !Land ocean mask at 1.0 deg (1:land, 0:ocean)
   real   ::ALT          !altitude (m above MSL)
   real   ::Albedo_soil0 !albedo, default
   real   ::W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real   ::W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real   ::W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real   ::W_mat        !matrix potential (m, -0.0001 -> -3.0)
   integer::SoilClass    !
   
!Others
   real    LAT, LON              !latitude and logitude for simulate
   integer GlobalZone            !ID number of global zone
   integer year, doy, dat, point !Counters
   integer i, j                  !for general usage
   
!_____________ Set Variables, for wide area computations
!For reading data
   real,dimension(Day_in_Year, YearMaxClimate1, 1:8)::dataREAD1  !1850~2005
   real,dimension(Day_in_Year, YearMaxClimate2, 1:8)::dataREAD2  !2006~2100
   
!For I/O
   character(len= 3) nam_lat
   character(len= 3) nam_lon
   character(len= 2) nam_dat
   integer i1, i2, i3
   integer file_no_grid1, file_no_grid2, file_no_grid3
   
!______________________Set varieties of reference number for this grid cell
!Device number for I/O
   file_no_grid1 = myid + 100 !For output files1, and for general usage
   file_no_grid2 = myid + 200 !For reading spinup files
   file_no_grid3 = myid + 300 !For writing spinup files
   
!Reference number for location
   ! LAT: north +, south - (decimalized)
   ! LON: east  +, west  - (decimalized)
   LAT   =    90.0 - (real(latNo)-0.5) * 0.25
   LON   = - 180.0 + (real(lonNo)-0.5) * 0.25
   point = (90-int(LAT)-1)*360 + int(LON+180) + 1 !grid point number @1.0 deg system
   
!Characters for designating latitude and longitude
   i1 = int(  latNo               /100 ) ; nam_lat(1:1) = char(i1+48)
   i2 = int( (latNo -i1*100)      /10  ) ; nam_lat(2:2) = char(i2+48)
   i3 =    (  latNo -i1*100-i2*10      ) ; nam_lat(3:3) = char(i3+48)
   
   i1 = int(  lonNo               /100 ) ; nam_lon(1:1) = char(i1+48)
   i2 = int( (lonNo -i1*100)      /10  ) ; nam_lon(2:2) = char(i2+48)
   i3 =    (  lonNo -i1*100-i2*10      ) ; nam_lon(3:3) = char(i3+48)
   
!GlobalZone: Set Location category
   if (LON>=-20 .and. 60>=LON .and. 23.0>=LAT ) then
      !African continent
      GlobalZone = 1
   elseif (LON>=100 .and. 170>=LON .and. 50.0<=LAT) then
      !Eastern Siberia
      GlobalZone = 2
   else
      !Default
      GlobalZone = 0
   endif
   
!______________ Read Location data
   open (file_no_grid1, file='land_prop.txt', status='OLD')
   do i=1, point
      read(file_no_grid1,*) Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt, SoilClass
   end do
   close(file_no_grid1)
   if (W_fi   > W_sat ) W_fi   = W_sat
   if (W_wilt > W_sat ) W_wilt = W_sat
   
!_____________ Prepare Climate Data
!Set sizes of allocatable climate data table
   allocate (tmp_air       (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !1
   allocate (tmp_air_range (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !2
   allocate (prec          (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !3
   allocate (rad_short     (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !4
   allocate (rad_long      (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !5
   allocate (wind          (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !6
   allocate (rh            (Day_in_Year, YearMaxClimate1+YearMaxClimate2)          ) !7
   allocate (tmp_soil      (Day_in_Year, YearMaxClimate1+YearMaxClimate2, NumSoil) ) ! 
   
   !Read climate data, 1850~2005
   i = 7 * 365 * YearMaxClimate1 * 4 !file size for each variables = 4byte
   Open (file_no_grid1, file=Loc_climate_data1//nam_lat//'/'//nam_lon//'.dat', &
         access='direct', recl=i, status='old')
      read (file_no_grid1,rec=1) &
      (((dataREAD1(doy,year,dat), dat=1,7), doy=1,365), year=1, YearMaxClimate1)
   Close(file_no_grid1)
   
   !Read climate data, 2006~2100
   i = 7 * 365 * YearMaxClimate2  * 4 !f[^TCY * 4byte
   !Open (file_no_grid1, file=Loc_climate_data2//nam_lat//'/'//nam_lon//'.dat', &
    !     access='direct', recl=i, status='old')
     ! read (file_no_grid1,rec=1) &
     ! (((dataREAD2(doy,year,dat), dat=1,7), doy=1,365), year=1, YearMaxClimate2)
  ! Close(file_no_grid1)
   
!_____________ Climate data Conversion
   do year=1,YearMaxClimate1
   do doy=1,Day_in_Year
      tmp_air      (doy,year)   = dataREAD1(doy,year,1)
      tmp_soil     (doy,year,:) = dataREAD1(doy,year,1) !substitute air-temp.
      tmp_air_range(doy,year)   = dataREAD1(doy,year,2)
      prec         (doy,year)   = dataREAD1(doy,year,3)
      rad_short    (doy,year)   = dataREAD1(doy,year,4)
      rad_long     (doy,year)   = dataREAD1(doy,year,5)
      wind         (doy,year)   = dataREAD1(doy,year,6)
      if (dataREAD1(doy,year,7) > 100) then
              rh   (doy,year) = 100
      else
              rh   (doy,year) = dataREAD1(doy,year,7)
      endif

   enddo
   enddo
   
   !do year = YearMaxClimate1+1, YearMaxClimate1+YearMaxClimate2
   !do doy=1,Day_in_Year
    !  tmp_air      (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,1)
     ! tmp_soil     (doy, year,:) = dataREAD2(doy,year-YearMaxClimate1,1) !substitute air-temp.
      !tmp_air_range(doy, year)   = dataREAD2(doy,year-YearMaxClimate1,2)
     ! prec         (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,3)
     ! rad_short    (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,4)
     ! rad_long     (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,5)
     ! wind         (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,6)
     ! rh           (doy, year)   = dataREAD2(doy,year-YearMaxClimate1,7)
  ! enddo
  ! enddo
   
!___________ For employing different random seed for each run (by Shigeki Ikeda @ Kyoto Univ.)
   IF (Flag_randomization) then
      call random_seed(size=seedsize)
      allocate(seed(seedsize))
      
      do size_count=1,seedsize
      call system_clock(count=clock)
      seed(size_count)=clock
      end do
      
      call random_seed(put=seed)
   EndIf
   
!_____________ !Call simulation loop for each grid
!Open I/O files
   !for writing output files
   if (Flag_output_write) then
      open (file_no_grid1, &
      file = Loc_result_files//nam_lat//'_'//nam_lon//'_out.txt')
   endif
   
   !for reading spinup files
   if (Flag_spinup_read) then
     open ( file_no_grid2, &
     file = Loc_spnin_files//nam_lat//'_'//nam_lon//'_spnin.dat', &
     access="sequential", form="unformatted", status="old", action="read")
   endif
   
   !for writing spinup files
   if (Flag_spinup_write) then
      open ( file_no_grid3, &
      file = Loc_result_files//nam_lat//'_'//nam_lon//'_spnout.dat', &
      access="sequential", form="unformatted", action="write")
   endif
   
!Call main program
   Call main_loop ( &
   LAT, LON, GlobalZone, YearMaxClimate1+YearMaxClimate2, YearMaxCO2, &
   file_no_grid1, file_no_grid2, file_no_grid3, &
   tmp_air(:,:), tmp_air_range(:,:), prec(:,:), rad_short(:,:), rad_long(:,:), &
   wind(:,:), rh(:,:), tmp_soil(:,:,:), &
   aco2_1850to2100(:), ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, SoilClass)
   
!Close I/O files
   if (Flag_output_write) close ( file_no_grid1 ) !æèÜÆßn}pÌoÍ
   if (Flag_spinup_read ) close ( file_no_grid2 ) !for reading spinup files
   if (Flag_spinup_write) close ( file_no_grid3 ) !for writing spinup files
   
END Subroutine start



!*************************************************************************************************
! Analysis output values of Last 10 yers
!*************************************************************************************************
Subroutine after_sim1(LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd, &
                      Loc_result_files, Loc_analysis_files, landmask)
   
   use data_structure, only: PFT_no
   implicit none
   
!_____________ Set Augment
   integer         ,intent(IN):: LatMax, LonMax, LatNoStart, LatNoEnd, LonNoStart, LonNoEnd
   character(len=*),intent(IN):: Loc_result_files, Loc_analysis_files
   
   logical,dimension(LonMax, LatMax),intent(IN):: landmask !Land Ocean mask
   
!_____________ Set parameters
   !Period for writing monthly outputs @ year
   !Monthly data is written for last YearForMean years of the simulation
   integer,parameter::YearForMean = 10 !This value has to be same as one in the main_mpi.f90
   
!_____________ Set variables
!Array for inputing result files
   real,dimension(1:LatMax, 1:LonMax):: Larch_area !Larch Area (fraction)
   
   integer,dimension(                                               2)::data1_read
   integer,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12,  2)::data1
   
   real   ,dimension(                                              28)::data2_read
   real   ,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12, 28)::data2

   real   ,dimension(                                              16)::data3_read
   real   ,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12, 16)::data3
 
   real,dimension(                                              16)::data4_read
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12, 16)::data4

   integer,dimension(                                               1)::data5_read
   integer,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 12,  1)::data5

!Variables for outputting geographic distribution
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd):: &
      out_water        , & !Water content @ top soil layer, annual average
      out_precipitation, & !Annual precip.
      out_gpp          , & !GPP (kg C/ m2/ year)
      out_npp          , & !NPP (kg C/ m2/ year)
      out_nep          , & !NEP (kg C/ m2/ year)
      out_csoil        , & !
      out_hr           , & !HR  (kg C/ m2/ year)
      out_lai_amean    , & !Annual mean of LAI (All   PFTs)
      out_lai_amean_t  , & !Annual mean of LAI (Woody PFTs)
      out_lai_amean_g  , & !Annual mean of LAI (Grass PFTs)
      out_lai_max      , & !Annula maximum of LAI (All   PFTs)
      out_lai_max_t    , & !Annula maximum of LAI (Woody PFTs)
      out_lai_max_g    , & !Annula maximum of LAI (Grass PFTs)
      out_ald_max      , & !Annula Maximum of Active Layer Depth (m)
      out_water_JJA        !Available soil water at top 5 layers during JJA(mm)
   
   real,dimension(LatNoStart:LatNoEnd, LonNoStart:LonNoEnd, 16)::&
      out_npp_pft      ,&  !NPP for PFT
      out_lai_pft          !lai for PFT
   !For summary statistics
   real sum_wbiomass  ! 1, Carbon in Woody biomass [kg C / m2]
   real sum_gbiomass  ! 2, Carbon in Grass biomass [kg C / m2]
   real sum_litter    ! 3, Carbon in litter        [kg C / m2]
   real sum_som_int   ! 4, Carbon in som_int       [kg C / m2]
   real sum_som_slow  ! 5, Carbon in som_slow      [kg C / m2]
   real sum_gpp       ! 9, GPP                     [kg C / m2 / month]
   real sum_npp       !10, NPP                     [kg C / m2 / month]
   real sum_nep       !11, NEP                     [kg C / m2 / month]
   real sum_hr        !12, Heterotrophic resp.     [kg C / m2 / month]
   real sum_runoff    !15, runoff                  [mm/month]
   real sum_intercept !16, interception            [mm/month]
   real sum_evapor    !17, evaporation             [mm/month]
   real sum_transpi   !18, transpiration           [mm/month]
   
   !etc
   character(len=7) fname            !For file name
   integer i, i1, i2, i3, i4, i5, i6    !For processing file name
   integer lat, lon, month, year, p     !Loop counter
   real    x, y, sum_weight
   real    a1, a2, a3, a4, a5, a6
   
!_______________ Read result files and make average value
write(*,*) "Reading result files and making output variables"
   !Initialize
   data1(:,:,:,:) =   0
   data2(:,:,:,:) = 0.0
   data3(:,:,:,:) = 0.0
   data4(:,:,:,:) =   0
   data5(:,:,:,:) =   0

   !Loop for each simulation grid
   DO lat = LatNoStart, LatNoEnd
   DO lon = LonNoStart, LonNoEnd
   IF (landmask(lon,lat) .eqv. .false.) cycle
      
      !Make letter string for designating a file
      i1 = int(  lat              /100 ) ; fname(1:1)=char(i1+48)
      i2 = int( (lat-i1*100)      / 10 ) ; fname(2:2)=char(i2+48)
      i3 =    (  lat-i1*100-i2*10      ) ; fname(3:3)=char(i3+48)
                                           fname(4:4)='_'
      i4 = int(  lon              /100 ) ; fname(5:5)=char(i4+48)
      i5 = int( (lon-i4*100)      / 10 ) ; fname(6:6)=char(i5+48)
      i6 =    (  lon-i4*100-i5*10      ) ; fname(7:7)=char(i6+48)
      
      !Read output files and make average values
      open (1, file=Loc_result_files//fname//'_out.txt', status='OLD')
         
         !Sumup all values for YearForMean year (except for Biome code)
         do year  =1, YearForMean
            do month =1, 12
               read(1,*) data1_read(:), data2_read(:), data3_read(:), data4_read(:), data5_read(:)
               data1(lat,lon,month,1) =                          data1_read(1) !Biome code
               data1(lat,lon,month,2) = data1(lat,lon,month,2) + data1_read(2) !Drought days
               data2(lat,lon,month,:) = data2(lat,lon,month,:) + data2_read(:) !Other Variables
               data3(lat,lon,month,:) = data3(lat,lon,month,:) + data3_read(:) !Other Variables
               data4(lat,lon,month,:) = data4(lat,lon,month,:) + data4_read(:) !Other Variables
               data5(lat,lon,month,1) =                          data5_read(1) !dominant PFT
            enddo
         enddo
         
         !Make average values by dividing YearForMean (except for Biome code)
         do month =1, 12
            data1(lat,lon,month,2) = data1(lat,lon,month,2) / real(YearForMean)
            data2(lat,lon,month,:) = data2(lat,lon,month,:) / real(YearForMean)
            data3(lat,lon,month,:) = data3(lat,lon,month,:) / real(YearForMean)
            data4(lat,lon,month,:) = data4(lat,lon,month,:) / real(YearForMean)
         enddo
         
      close (1)
      
   END DO
   END DO
   
!_______________ Prepare output variables
write(*,*) "Prepare output variables"
   !initialize
   out_water         (:,:) = 0.0
   out_precipitation (:,:) = 0.0
   out_csoil         (:,:) = 0.0
   
   out_gpp           (:,:) = 0.0
   out_npp           (:,:) = 0.0
   out_nep           (:,:) = 0.0
   out_hr            (:,:) = 0.0
   
   out_lai_amean     (:,:) = 0.0
   out_lai_amean_t   (:,:) = 0.0
   out_lai_amean_g   (:,:) = 0.0
   
   out_lai_max       (:,:) = 0.0
   out_lai_max_t     (:,:) = 0.0
   out_lai_max_g     (:,:) = 0.0
   
   out_ald_max       (:,:) = 0.0
   out_water_JJA     (:,:) = 0.0
   out_npp_pft       (:,:,:) = 0.0   
   out_lai_pft       (:,:,:) = 0.0

   Do lat   = LatNoStart, LatNoEnd
   Do lon   = LonNoStart, LonNoEnd
      
      Do month = 1, 12
      !Soil water contents
      out_water(lat,lon)= out_water(lat,lon)+ data2(lat,lon,month, 6) / 12.0
      
      !Annual GPP
      out_gpp  (lat,lon)= out_gpp  (lat,lon)+ data2(lat,lon,month, 9)
      
      !Annual NPP
      out_npp  (lat,lon)= out_npp  (lat,lon)+ data2(lat,lon,month,10)
      
      !Annual NEP
      out_nep  (lat,lon)= out_nep  (lat,lon)+ data2(lat,lon,month,11)
      
      !Annual Heterotrophic resp.
      out_hr   (lat,lon)= out_hr   (lat,lon)+ data2(lat,lon,month,25)
      
      !Annual means of LAI
      out_lai_amean  (lat,lon) = out_lai_amean  (lat,lon) + data2(lat,lon,month,13) / 12.0 &
                                                          + data2(lat,lon,month,14) / 12.0
      out_lai_amean_t(lat,lon) = out_lai_amean_t(lat,lon) + data2(lat,lon,month,13) / 12.0
      out_lai_amean_g(lat,lon) = out_lai_amean_g(lat,lon) + data2(lat,lon,month,14) / 12.0
      
      !Annual maximums of LAI
      out_lai_max  (lat,lon) = max( out_lai_max  (lat,lon), data2(lat,lon,month,13) &
                                                          + data2(lat,lon,month,14) )
      out_lai_max_t(lat,lon) = max( out_lai_max_t(lat,lon), data2(lat,lon,month,13) )
      out_lai_max_g(lat,lon) = max( out_lai_max_g(lat,lon), data2(lat,lon,month,14) )
      
      !Annual Precipitation
      out_precipitation(lat,lon) = out_precipitation(lat,lon)+ sum(data2(lat,lon,month,15:18))
      
      !Hydorology Related
      out_ald_max   (lat,lon) = max ( out_ald_max(lat,lon), data2(lat,lon,month,27) )
      
      End Do
      
      !Hydorology Related
      Do month = 6, 8
      out_water_JJA (lat,lon) = out_water_JJA (lat,lon) + data2(lat,lon,month,28) / 3.0
      End Do
     
      do p=1, PFT_no 
        !NPP for PFT
        out_npp_pft (lat,lon,p) = out_npp_pft (lat,lon,p) + data3(lat,lon,month,p) ! sumNPPPFT is already the yearly sum of npp in etc.f90
      enddo

      !lai for PFT
      !Do month = 1, 12
      out_lai_pft (lat,lon,:) = data4(lat,lon,12,:)
      !End Do
   End Do
   End Do
   
!_______________ Write output files 1
write(*,*) "Writing result maps 1"
   Open ( 1, file=Loc_analysis_files//'out_biome.txt'        )
   Open ( 2, file=Loc_analysis_files//'out_wbiomass.txt'     )
   Open ( 3, file=Loc_analysis_files//'out_water1.txt'       )
   Open ( 4, file=Loc_analysis_files//'out_gpp.txt'          )
   Open ( 5, file=Loc_analysis_files//'out_npp.txt'          )
   Open ( 6, file=Loc_analysis_files//'out_nep.txt'          )
   Open ( 7, file=Loc_analysis_files//'out_hr.txt'           )
   
   Open ( 8, file=Loc_analysis_files//'out_precipitation.txt')
   Open ( 9, file=Loc_analysis_files//'out_fire.txt'         )
   
   Open (10, file=Loc_analysis_files//'out_lai_max.txt'      )
   Open (11, file=Loc_analysis_files//'out_lai_max_t.txt'    )
   Open (12, file=Loc_analysis_files//'out_lai_max_g.txt'    )
   
   Open (20, file=Loc_analysis_files//'out_lai_amean.txt'    )
   Open (21, file=Loc_analysis_files//'out_lai_amean_t.txt'  )
   Open (22, file=Loc_analysis_files//'out_lai_amean_g.txt'  )
   
   Open (23, file=Loc_analysis_files//'out_ald_max.txt'      )
   Open (24, file=Loc_analysis_files//'out_water_JJA.txt'    )
   
   Open (25, file=Loc_analysis_files//'out_csoil.txt'        )
   Open (26, file=Loc_analysis_files//'out_npppft01.txt'     )
   Open (27, file=Loc_analysis_files//'out_npppft02.txt'     )
   Open (28, file=Loc_analysis_files//'out_npppft03.txt'     )
   Open (29, file=Loc_analysis_files//'out_npppft04.txt'     )
   Open (30, file=Loc_analysis_files//'out_npppft05.txt'     )
   Open (31, file=Loc_analysis_files//'out_npppft06.txt'     )
   Open (32, file=Loc_analysis_files//'out_npppft07.txt'     )
   Open (33, file=Loc_analysis_files//'out_npppft08.txt'     )
   Open (34, file=Loc_analysis_files//'out_npppft09.txt'     )
   Open (35, file=Loc_analysis_files//'out_npppft10.txt'     )
   Open (36, file=Loc_analysis_files//'out_npppft11.txt'     )
   Open (37, file=Loc_analysis_files//'out_npppft12.txt'     )
   Open (38, file=Loc_analysis_files//'out_npppft13.txt'     )
   Open (39, file=Loc_analysis_files//'out_npppft14.txt'     )
   Open (40, file=Loc_analysis_files//'out_npppft15.txt'     )
   Open (41, file=Loc_analysis_files//'out_npppft16.txt'     )
   
   Open (42, file=Loc_analysis_files//'out_laipft01.txt'     )
   Open (43, file=Loc_analysis_files//'out_laipft02.txt'     )
   Open (44, file=Loc_analysis_files//'out_laipft03.txt'     )
   Open (45, file=Loc_analysis_files//'out_laipft04.txt'     )
   Open (46, file=Loc_analysis_files//'out_laipft05.txt'     )
   Open (47, file=Loc_analysis_files//'out_laipft06.txt'     )
   Open (48, file=Loc_analysis_files//'out_laipft07.txt'     )
   Open (49, file=Loc_analysis_files//'out_laipft08.txt'     )
   Open (50, file=Loc_analysis_files//'out_laipft09.txt'     )
   Open (51, file=Loc_analysis_files//'out_laipft10.txt'     )
   Open (52, file=Loc_analysis_files//'out_laipft11.txt'     )
   Open (53, file=Loc_analysis_files//'out_laipft12.txt'     )
   Open (54, file=Loc_analysis_files//'out_laipft13.txt'     )
   Open (55, file=Loc_analysis_files//'out_laipft14.txt'     )
   Open (56, file=Loc_analysis_files//'out_laipft15.txt'     )
   Open (57, file=Loc_analysis_files//'out_laipft16.txt'     )
   Open (58, file=Loc_analysis_files//'out_pftdominant.txt'  )
   Do lat = LatNoStart, LatNoEnd
      Do lon = LonNoStart, LonNoEnd
      write( 1,'(  i2,a)', advance='no') data1     (lat,lon,12, 1),  ',' !Biome
      write( 2,'(f7.2,a)', advance='no') data2     (lat,lon,12, 1),  ',' !Woody biomass
      write( 3,'(f9.1,a)', advance='no') out_water (lat,lon)      ,  ',' !Water content @ top soil layer, annual average
      write( 4,'(f8.3,a)', advance='no') out_gpp          (lat,lon), ',' !GPP (kg C/ m2/ year)
      write( 5,'(f8.3,a)', advance='no') out_npp          (lat,lon), ',' !NPP (kg C/ m2/ year)
      write( 6,'(f8.3,a)', advance='no') out_nep          (lat,lon), ',' !NEP (kg C/ m2/ year)
      write( 7,'(f8.3,a)', advance='no') out_hr           (lat,lon), ',' !HR  (kg C/ m2/ year)
      
      write( 8,'(f9.1,a)', advance='no') out_precipitation(lat,lon      ), ',' !Annual precip.
      write( 9,'(f5.3,a)', advance='no') data2            (lat,lon,12,21), ',' !Fire frequency
      
      !Annual maximum LAI
      write(10,'(f4.1,a)', advance='no') out_lai_max      (lat,lon), ',' !All PFTs
      write(11,'(f4.1,a)', advance='no') out_lai_max_t    (lat,lon), ',' !Woody PFTs
      write(12,'(f4.1,a)', advance='no') out_lai_max_g    (lat,lon), ',' !Grass PFTs
      
      !Annual average LAI
      write(20,'(f4.1,a)', advance='no') out_lai_amean    (lat,lon), ',' !All PFTs
      write(21,'(f4.1,a)', advance='no') out_lai_amean_t  (lat,lon), ',' !Woody PFTs
      write(22,'(f4.1,a)', advance='no') out_lai_amean_g  (lat,lon), ',' !Grass PFTs
      
      !Annual maximum of activeaverage LAI
      write(23,'(f5.3,a)', advance='no') out_ald_max      (lat,lon), ',' !Active Layer Maximum
      
      !Annual maximum of activeaverage LAI
      write(24,'(f5.1,a)', advance='no') out_water_JJA    (lat,lon), ',' !JJA available water [mm]
      
      !Soil Carbon
      write(25,'(f9.2,a)', advance='no') sum(data2(lat,lon,12,3:5)),  ',' !Soil Carbon [kg C / m2]
 
      do p=0,PFT_no-1

      !NPP for PFTs
      write(26+p,'(f12.3,a)', advance='no') out_npp_pft(lat,lon,p+1),  ',' !npp for pft
      
      !lai for PFTs
      write(42+p,'(f12.3,a)', advance='no') out_lai_pft(lat,lon,p+1),  ',' !lai for pft
      enddo

      !dominant pft
      write( 58,'(  i2,a)', advance='no') data5     (lat,lon,12, 1),  ',' !dominant pft

      End Do
      
      !Insert feed code
      write( 1,*); write( 2,*); write( 3,*); write( 4,*); write( 5,*); write( 6,*); write( 7,*)
      write( 8,*); write( 9,*)
      write(10,*); write(11,*); write(12,*)
      write(20,*); write(21,*); write(22,*)
      write(23,*); write(24,*); write(25,*)
      Do p =26,58
      write(p,*)
      Enddo 
      
   End Do
   
   Close ( 1); Close ( 2); Close ( 3); Close ( 4); Close ( 5); Close ( 6); Close ( 7)
   Close ( 8); Close ( 9)
   Close (10); Close (11); Close (12)
   Close (20); Close (21); Close (22)
   Close (23); Close (24); Close (25)
   Do p =26,58
   Close (p)
   Enddo 
END Subroutine after_sim1



