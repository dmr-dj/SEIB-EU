#*********************************************************************
# R code for drawing map output variables of SEIB-DGVM ver. 3.0-      
#                                                                     
# All rights are reserved by Dr. Hisashi SATO @ JAMSTEC               
#*********************************************************************

#____________________ Common procedure ____________________ 
#Set environment
   #Working directory
   setwd('./../result_visualize/') 
   
   #Missing value
   Missing  =  0
   
   #Grid size (length on a side @ deg)
   GridSize   = 0.25
   
   #Pixel size of a grid
   PixSize   = 2.0
   
   #Visualizing area (designated by grid numbers)
   LatNoStart = 73
   LatNoEnd   = 228
   LonNoStart = 665
   LonNoEnd   = 920
   
#Compute valiables for coodination and image size
   
   #Row and column numbers of drawing area
   Lat      = LatNoEnd - LatNoStart + 1
   Lon      = LonNoEnd - LonNoStart + 1
   
   #Vertical and horizontal picture size @ pixel
   width_size  = Lon * PixSize
   height_size = Lat * PixSize
   
   #Latitude and longitude at the center of grid of the most south-west grid cell
   Lat1_loc =    90.0 - (LatNoEnd   - 0.5) * GridSize #(North:+ , Sourth:-)
   Lon1_loc = - 180.0 + (LonNoStart - 0.5) * GridSize #(West :- , East  :+)
   
   #Latitude coordinate
   y <- array(0.0, dim=c(Lat))                          #Prepare array y
   for(i in 1:Lat)   y[i] <- Lat1_loc + (i-1)*GridSize  #Input latitude from south to north
   
   #Longitude coordinate
   x <- array(0.0, dim=c(Lon))                          #Prepare attay x
   for(i in 1:Lon)   x[i] <- Lon1_loc + (i-1)*GridSize  #Input longitude from west to east
   
#Other preparation
   #Activate map library
   library(maps)
   library(RColorBrewer)
  #library(colorRamps)
   
   #Prepare data aray
   z  <- array(Missing, dim=c(Lon,Lat))

#____________________ Subroutines for color palette ____________________ 
set_color_topo <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(topo.colors(num_col))    #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "white"                 #Change minimum class color
   return(col)                       }
   
set_color_topo2 <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(topo.colors(num_col))    #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]        = "white"           #Change minimum class color
   col[num_col]  = "gray "           #Change maximum class color
   return(col)                       }
   
set_color_heat <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(heat.colors(num_col))    #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "white"                 #Change minimum class color
   return(col)                       }
   
set_color_terrain <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(terrain.colors(num_col)) #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "antiquewhite1"         #Change minimum class color
   return(col)                       }

set_color_PlusMinus <- function(num_col) {
   col <-  c(num_col)                  #Prepare color variable (this must be division number - 1)
#  col <- c(matlab.like(num_col))      #Give color from a palette
#  col <- c(cm.colors(num_col))        #Give color from a palette
   col <- brewer.pal(num_col, "PiYG")  #Give color from a palette
   col[6] <- "white"                   #For Zero Color
#  col <- col[length(col):1]           #Turn upside down for color palette
   return(col)                       }


#____________________ Subroutines for data reading and coordination conversion ____________________ 
read_data <- function(fname) {
   d <- read.csv(fname, header=F)
   for (i in 1:Lat) {
   for (j in 1:Lon) {
      #Tuen upside down for latitude, and replace row and column
      z[j,i] <- d[Lat-i+1,j] 
   }
   }
   return(z)
}

#____________________ Subroutine for drawing color pannel ____________________
# LabelName :Label on the top of color pannel
# DivedNum  :Cell number of color pannel
# LabelNum  :Number of value below the color pannel
# IncreStep :Increment of color pannel

draw_panel <- function(LabelName, DivedNum, LabelNum, IncreStep) {
   x_start <-   -13
   x_width <-    2
   y_start <-   73 #�]�����܂�
   y_width <-    2
   
   # �g�O�ւ̕`�������
   par(xpd=T)
   
   #Write label on top of the color pannel
   text(x_start+9, y_start-2*y_width, pos=1, LabelName)
   
   #Draw white belt under color pannel
   x_end <-  x_start+x_width*LabelNum
   polygon( c(x_start, x_end, x_end, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col='white') 
   
   #Draw color pannel
   if (DivedNum==0 || LabelNum==0) {return}
   i <- floor(DivedNum/ LabelNum)
   for (j in 1:DivedNum) {
      polygon( c(x_start, x_start+x_width, x_start+x_width, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col=col[j])
      if (floor(j/i) == j/i) {text(x_start+0.5*x_width, y_start-y_width, pos=1, j*IncreStep)}
      x_start <- x_start+x_width
   
   }
}


#____________________ Subroutine for drawing a distribution map 1 ____________________
draw_dist <- function(DivedNum,PannelStep,col,x,y,z) {
   PannelMax  <- DivedNum * PannelStep                #�J���[�p�l���̍ő�l���Z�o
   br  <- c(seq(from=0, to=PannelMax, by=PannelStep)) #�l�̕����_�̐ݒ�P(�����̎d�؂�)
   br[DivedNum+1]  = PannelMax*10                     #�l�̕����_�̐ݒ�Q(�V��)
   par(mar = c(0,0,0,0) )                             #���E���E��E�E�̏��œ����}�[�W����ݒ�
   frame()                                            #��ʂ̃N���A
   image(x, y, z, breaks=br, col=col, xlab='', ylab='', axes = FALSE) #�}�b�s���O
   map(add=T, interior=FALSE)                                         #�n�}���d�ˏ���
   }
   
#____________________ Subroutine for drawing a distribution map 2 ____________________
draw_dist_LAI <- function(DivedNum,PannelStep,col,x,y,z) {
   PannelMax  <- DivedNum * PannelStep                #�J���[�p�l���̍ő�l���Z�o
   br  <- c(seq(from=0, to=PannelMax, by=PannelStep)) #�l�̕����_�̐ݒ�P(�����̎d�؂�)
   br[DivedNum+1]  = PannelMax*10                     #�l�̕����_�̐ݒ�Q(�V��)
   br[1]           = 0.01                             #�l�̕����_�̐ݒ�R(�n�_)
   par(mar = c(0,0,0,0) )                             #���E���E��E�E�̏��œ����}�[�W����ݒ�
   frame()                                            #��ʂ̃N���A
   image(x, y, z, breaks=br, col=col, xlab='', ylab='', axes = FALSE) #�}�b�s���O
   map(add=T, interior=FALSE)                                         #�n�}���d�ˏ���
   }
   
#____________________ Biome ____________________ 
   DivedNum   <- 12 #�J���[�p�l���̕�����
   PannelStep <- 1  #�J���[�p�l���̑���
   
   z <- read_data('out_biome.txt') #�f�[�^�ǂ݂����ƁA���`
   
   #�F�̐ݒ���s��
   col <- set_color_topo(DivedNum) 
   col[ 0] <- "white"         # 0: - water -
   col[ 1] <- "lightblue"         # 1: Polar desert
   col[ 2] <- "plum1"  # 2: Arctic/Alpine-tundra
   col[ 3] <- "red"           # 3: Tropical evergreen forest
   col[ 4] <- "purple"        # 4: Tropical deciduous forest
   col[ 5] <- "yellowgreen"   # 5: Temperate conifer forest
   col[ 6] <- "yellow4"       # 6: Temperate broad-leaved evergreen forest
   col[ 7] <- "limegreen"     # 7: Temperate deciduous forest
   col[ 8] <- "aquamarine4"   # 8: Boreal evergreen forest
   col[ 9] <- "steelblue1"    # 9: Boreal deciduous forest
   col[10] <- "steelblue4"    #10: Xeric woodland / scrub
   col[11] <- "bisque"        #11: Grassland / steppe / Savanna
   col[12] <- "orange"       #12: Desert
   
   png('out_biome.png', width=width_size, height=height_size) #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   
   #�}��̕`�� ______
   x_start <-   -13
   x_width <-    2
   y_start <-    73  #�]�����܂�
   y_width <-    2
   
   # �g�O�ւ̕`�������
   par(xpd=T)
   
   #Write label on top of the color pannel
   text(x_start+5, y_start-2*y_width, pos=1, "Biome Type")
   
   #Draw color pannel
   for (j in 1:DivedNum) {
      polygon( c(x_start, x_start+x_width, x_start+x_width, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col=col[j])
      text(x_start+0.5*x_width, y_start-y_width, pos=1, j)
      x_start <- x_start+x_width
   }
   
   #�f�o�C�X�h���C�o���� ______
   dev.off()
   
#____________________ Fire number (n/year) ____________________ 
   DivedNum   <- 20      #�J���[�p�l���̕�����
   PannelStep <- 0.0025  #�J���[�p�l���̑���
   
   z   <- read_data('out_fire.txt')                               #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum)                                #�F�̐ݒ���s��
   
   png('out_fire.png', width=width_size, height=height_size)      #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)                 #�}�{�̂̕`��
   draw_panel('Fire Frequency (n/year)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                      #�f�o�C�X�h���C�o����
   
#____________________ Biomass ____________________ 
   DivedNum   <- 20   #�J���[�p�l���̕�����
   PannelStep <- 1.0 #�J���[�p�l���̑���
   
   z <- read_data('out_wbiomass.txt')                            #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_heat(DivedNum)                               #�F�̐ݒ���s��
   
   png('out_wbiomass.png', width=width_size, height=height_size) #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)                #�}�{�̂̕`��
   draw_panel('Biomass (KgC/m2)', DivedNum, 4, PannelStep)       #�J���[�p�l�����d�ˏ���
   dev.off()                                                     #�f�o�C�X�h���C�o����

#____________________ GPP ____________________ 
   DivedNum   <- 12                    #�J���[�p�l���̕�����
   PannelStep <- 0.2                   #�J���[�p�l���̑���
   
   z <- read_data('out_gpp.txt')       #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum)     #�F�̐ݒ���s��
   
   png('out_gpp.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('GPP (gC/m2/yr)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off()                                                  #�f�o�C�X�h���C�o����
   
#____________________ NPP ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npp.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npp.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off()                                                  #�f�o�C�X�h���C�o����
   
#____________________ NEP ____________________ 
   DivedNum   <- 11               #�J���[�p�l���̕����� (��ɂȂ�悤�ɐݒ肷�邱��)
   PannelStep <- 0.1              #�J���[�p�l���̑���
   
   z <- read_data('out_nep.txt')         #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_PlusMinus(DivedNum)  #�F�̐ݒ���s��
   
   png('out_nep.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   par(xpd=T)                                                 #�g�O�ւ̕`�������
   
   PannelMin  <- -0.5*(DivedNum-1) * PannelStep - 0.5 * PannelStep #�J���[�p�l���̍ő�l���Z�o
   PannelMax  <-  0.5*(DivedNum-1) * PannelStep + 0.5 * PannelStep #�J���[�p�l���̍ő�l���Z�o
   
   br  <- c(seq(from=PannelMin, to=PannelMax, by=PannelStep)) #�l�̕����_�̐ݒ�P(�����̎d�؂�)
   br[1]           = -PannelMax*10                            #�l�̕����_�̐ݒ�Q(�{�g��)
   br[DivedNum+1]  =  PannelMax*10                            #�l�̕����_�̐ݒ�R(�V��)
   
   par(mar = c(0,0,0.0,0) )                                           #���E���E��E�E�̏��œ����}�[�W����ݒ�
   frame()                                                            #��ʂ̃N���A
   image(x, y, z, breaks=br, col=col, xlab='', ylab='', axes = FALSE) #�}�b�s���O
   map(add=T, interior=FALSE)                                         #�n�}���d�ˏ���
   
#   #�J���[�p�l�����d�ˏ���
#   x_start <-   100
#   y_start <-   86 #�]�������ւ̂͂ݏo�����l�����Ă���
#   x_width <-    6
#   y_width <-    4
#   
#   LabelNum <-   5
#   IncreStep <-  PannelStep*1000
#   
#   text(x_start, y_start+2, pos=4, 'NEP (gC/m2/yr)') #Write label on top of the color pannel
#   
#   #Draw white belt under color pannel
#   x_end <-  x_start+x_width*DivedNum
#   polygon( c(x_start, x_end, x_end, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col='white')
#   
#   #Draw color pannel
#   i <- floor(DivedNum/ LabelNum)
#   for (j in 1:DivedNum) {
#      polygon( c(x_start, x_start+x_width, x_start+x_width, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col=col[j])
#      if (floor(j/i) == j/i) {text(x_start+0.5*x_width, y_start-y_width*0.7, pos=1, j*IncreStep-600)}
#      x_start <- x_start+x_width
#   }
   
   #�f�o�C�X�h���C�o����
   dev.off()                                                  
   
#____________________ HR ____________________ 
   DivedNum   <- 12               #�J���[�p�l���̕�����
   PannelStep <- 0.1              #�J���[�p�l���̑���
   
   z <- read_data('out_hr.txt')    #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_hr.png', width=width_size, height=height_size)    #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('HR (gC/m2/yr)', DivedNum, 6, PannelStep*1000)  #�J���[�p�l�����d�ˏ���
   dev.off()                                                  #�f�o�C�X�h���C�o����
   
#____________________ Maximum ALD (m) ____________________ 
   DivedNum   <- 20                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_ald_max.txt') #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo2(DivedNum)   #�F�̐ݒ���s��
   
   png('out_ald_max.png', width=width_size, height=height_size)           #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)                         #�}�{�̂̕`��
   draw_panel('Maximuum Active Layer Depth (m)', DivedNum, 5, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                              #�f�o�C�X�h���C�o����
   
#____________________ JJA Available Water within the Soil layers 1 to 5 ____________________ 
   DivedNum   <- 20                 #�J���[�p�l���̕�����
   PannelStep <- 0.02               #�J���[�p�l���̑���
   
   z <- read_data('out_water_JJA.txt') #�f�[�^�ǂ݂����ƁA���`
   for (i in 1:Lat) {
   for (j in 1:Lon) {
      z[j,i] <- z[j,i] / 500.0   #mm -> fraction
   }
   }
   col <- set_color_topo(DivedNum)  #�F�̐ݒ���s��
   
   png('out_water_JJA.png', width=width_size, height=height_size)                      #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)                                      #�}�{�̂̕`��
   draw_panel('JJA available Soil Water @ 0-50cm depth (fraction)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                                           #�f�o�C�X�h���C�o����
   
#____________________ LAI ____________________ 
   DivedNum   <- 14                   #�J���[�p�l���̕�����
   PannelStep <- 0.5                  #�J���[�p�l���̑���
   
   col <- set_color_terrain(DivedNum)  #�����_�ƐF�̐ݒ�
   
   #�`��a1
   z <- read_data('out_lai_amean.txt')                             #�f�[�^�ǂ݂����ƁA���`
   png('out_lai_amean.png', width=width_size, height=height_size)  #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #�}�{�̂̕`��
   draw_panel('Annual Mean LAI (m2/m2)', DivedNum, 7, PannelStep)  #�J���[�p�l�����d�ˏ���
   dev.off()                                                       #�f�o�C�X�h���C�o����
   
   #�`��a2
   z <- read_data('out_lai_amean_t.txt')                                        #�f�[�^�ǂ݂����ƁA���`
   png('out_lai_amean_t.png', width=width_size, height=height_size)             #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                           #�}�{�̂̕`��
   draw_panel('Annual Mean LAI of Woody PFTs(m2/m2)', DivedNum, 7, PannelStep)  #�J���[�p�l�����d�ˏ���
   dev.off()                                                                    #�f�o�C�X�h���C�o����
   
   #�`��a3
   z <- read_data('out_lai_amean_g.txt')                                        #�f�[�^�ǂ݂����ƁA���`
   png('out_lai_amean_g.png', width=width_size, height=height_size)             #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                           #�}�{�̂̕`��
   draw_panel('Annual Mean LAI of Grass PFT (m2/m2)', DivedNum, 7, PannelStep)  #�J���[�p�l�����d�ˏ���
   dev.off()                                                                    #�f�o�C�X�h���C�o����
   
   #�`��b1
   z <- read_data('out_lai_max.txt')                                 #�f�[�^�ǂ݂����ƁA���`
   png('out_lai_max.png', width=width_size, height=height_size)      #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                #�}�{�̂̕`��
   draw_panel('Annual Maximum LAI (m2/m2)', DivedNum, 7, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                         #�f�o�C�X�h���C�o����
   
   #�`��b2
   z <- read_data('out_lai_max_t.txt')                                            #�f�[�^�ǂ݂����ƁA���`
   png('out_lai_max_t.png', width=width_size, height=height_size)                 #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                             #�}�{�̂̕`��
   draw_panel('Annual Maximum LAI of tree PFTs (m2/m2)', DivedNum, 7, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                                      #�f�o�C�X�h���C�o����
   
   #�`��b3
   z  <- read_data('out_lai_max_g.txt')
   png('out_lai_max_g.png', width=width_size, height=height_size)                  #�f�o�C�X�h���C�o�J��
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                              #�}�{�̂̕`��
   draw_panel('Annual Maximum LAI of grass PFTs (m2/m2)', DivedNum, 7, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                                       #�f�o�C�X�h���C�o����
 #____________________ NPP_PFT01 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft01.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft01.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved evergreen     (TrBE_1)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off()                                                  #�f�o�C�X�h���C�o���� 
 #____________________ NPP_PFT02 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft02.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft02.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved evergreen     (TrBE_2)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT03 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft03.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft03.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved evergreen     (TrBE_3)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT04 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft04.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft04.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved evergreen     (TrBE_4)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT05 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft05.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft05.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved evergreen     (TrBE_5_Africa)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT06 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft06.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft06.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical broad-leaved raingreen     (TrBR_Africa)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT07 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft07.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft07.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Temperate needle-leaved evergreen   (TeNE)', DivedNum, 4, PannelStep*1) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT08 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft08.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft08.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Temperate broad-leaved evergreen    (TeBE)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT09 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft09.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft09.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Temperate broad-leaved summergreen  (TeBS)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT10 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft10.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft10.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Boreal needle-leaved evergreen      (BoNE)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT11 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft11.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft11.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Pice obovata (East Siberia)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT12 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft12.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft12.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Pinus sylvestris (East Siberia)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT13 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft13.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft13.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Boreal needle-leaved summergreen    (BoNS)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT14 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft14.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft14.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Boreal broad-leaved summergreen     (BoBS)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT15 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft15.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft15.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Temperate herbaceous(C3)            (TeH)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT16 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_npppft16.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_npppft16.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('NPP (gC/m2/yr) for Tropical herbaceous(C4)             (TrH)', DivedNum, 4, PannelStep*1000) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT01 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft01.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft01.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved evergreen     (TrBE_1)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off()                                                  #�f�o�C�X�h���C�o���� 
 #____________________ LAI_PFT02 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft02.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft02.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved evergreen     (TrBE_2)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT03 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft03.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft03.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved evergreen     (TrBE_3)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT04 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft04.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft04.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved evergreen     (TrBE_4)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT05 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft05.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft05.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved evergreen     (TrBE_5_Africa)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT06 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft06.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft06.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical broad-leaved raingreen     (TrBR_Africa)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT07 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.5               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft07.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft07.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Temperate needle-leaved evergreen   (TeNE)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT08 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.5               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft08.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft08.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Temperate broad-leaved evergreen    (TeBE)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT09 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.5               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft09.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft09.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Temperate broad-leaved summergreen  (TeBS)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT10 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.5               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft10.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft10.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Boreal needle-leaved evergreen      (BoNE)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT11 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft11.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft11.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Pice obovata (East Siberia)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT12 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft12.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft12.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Pinus sylvestris (East Siberia)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT13 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft13.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft13.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Boreal needle-leaved summergreen    (BoNS)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT14 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft14.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft14.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Boreal broad-leaved summergreen     (BoBS)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ LAI_PFT15 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft15.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft15.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Temperate herbaceous(C3)            (TeH)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 
 #____________________ NPP_PFT16 ____________________ 
   DivedNum   <- 12                #�J���[�p�l���̕�����
   PannelStep <- 0.1               #�J���[�p�l���̑���
   
   z <- read_data('out_laipft16.txt')   #�f�[�^�ǂ݂����ƁA���`
   col <- set_color_topo(DivedNum) #�F�̐ݒ���s��
   
   png('out_laipft16.png', width=width_size, height=height_size)   #�f�o�C�X�h���C�o�J��
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #�}�{�̂̕`��
   draw_panel('LAI max on a day for the year [m2/m2] for Tropical herbaceous(C4)             (TrH)', DivedNum, 4, PannelStep) #�J���[�p�l�����d�ˏ���
   dev.off() 