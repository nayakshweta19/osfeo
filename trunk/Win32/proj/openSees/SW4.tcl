# Create ModelBuilder with 2 dimensions and 2 DOF/node
wipe;
model basic -ndm 2 -ndf 2

# define UNITS
set in 1.; # define basic units — output units
set kip 1.; # define basic units — output units
set sec 1.; # define basic units - output units
set ft  [expr 12.*$in]; # define engineering units
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set lbf [expr $psi*$in*$in]; # pounds force
set pcf [expr $lbf/pow($ft,3)]; # pounds per cubic foot
set in2 [expr $in*$in]; # inchA2
set in4 [expr $in*$in*$in*$in]; # inchA4
set cm [expr $in/2.54]; # centimeter, needed for displacement input in
set PI [expr 2*asin(1.0)]; # define constants
set g [expr 32.2*$ft/pow($sec,2)]; # gravitational acceleration

# define input displacement history
set DHfile1 "TimeSeriesY.thf";

## create the material
## MATERIAL parameters
# nominal concrete compressive strength

set CompressiveStress [expr -5.35*$ksi];
set fc [expr $CompressiveStress];
set Ec [expr 57*sqrt(-$fc*1000)];

# unconfmed concrete
set fclU $fc; # UNCONFINED concrete (todeschini parabolic model), maximum stress
set epslU -0.003; # strain at maximum strength of unconfmed concrete
set fc2U [expr 0.2*$fclU]; # ultimate stress
set eps2U -0.01; # strain at ultimate stress
set lambda 0.1; # ratio between unloading slope at $eps2 and initial slope $Ec

# tensile-strength properties
set ftU [expr -0.14*$fclU]; # tensile strength +tension
set Ets [expr $ftU/0.002]; # tension softening stiffness

# DEFINE MATERIAL tags
set IDMat 1;
set IDMatFlg 2;
set IDMatSlab 3;

# back-bone steel stress-strain curve parameters for all materials
set Fy [expr 60*$ksi]; # STEEL yield stress
set Es [expr 29000*$ksi]; # modulus of steel
set epsY [expr $Fy/$Es]; # steel yield strain
set epsSH [expr 8*$epsY];
set Esh [expr 0.02*$Es]; # tangent stiffness at onset of StrainHardening
set Bs 0.01; # strain-hardening ratio
set Fy1 [expr 1.01*$Fy]; # steel stress post-yield
set epsY1 [expr $epsY+($Fy1-$Fy)/($Bs*$Es)]; # steel strain post-yield
set Fu [expr 1.16*$Fy];
set epsU 1.0; # ultimate strain of steel#

# hysteretic steel material parameters -- baseline
set pinchX 1.0; # pinching parameter for hysteretic model
set pinchY 1.0; # pinching parameter for hysteretic model

set damage1 0.0; # damage parameter for hysteretic model
set damage2 0.0; # damage parameter for hysteretic model
set betaMUsteel 0.0; # degraded unloading stiffness for hysteretic material based on MUA(-beta)

# steel02 and steel03 parameters — baseline
set R0 18;  # control the transition from elastic to plastic branches,
set cR1 0.925; # control the transition from elastic to plastic branches,
set cR2 0.15; # control the transition from elastic to plastic branches
set a2 0.1; # isotropic hardening parameter, associated with al
set a1 [expr $a2*($Fy/$Es)]; # isotropic hardening parameter, increase of comp.
                            # yield envelope as proportion of yield strength after a plastic strain
set a4 0.1; #ft isotropic hardening parameter, associated with a3
set a3 [expr $a4*($Fy/$Es)]; # isotropic hardening parameter, increase of tension yield envelope
                             # as proportion of yield strength after a plastic strain

# Reinforcing direction, D, and ratio, R, and diameter, d.
set rD1 0.0;
set rR1 0.0039;
set rdmm1 0.2362; #in

set rD2 90.0;
set rR2 0.0050;
set rdmm2 0.2362; #in

set FrR1 0.0123;
set FrR2 0.0686;
set Frd1 0.2362;
set Frd2 0.4724;

set SR1 0.003;
set SR2 0.06;

#MATERIAL assignments

nDMaterial NonlinearBS $IDMat $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $rR1 $rdmm1 $rD2 $rR2 $rdmm2 $Fy \
						$epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

nDMaterial NonlinearBS $IDMatFlg $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $FrR1 $Frd1 $rD2 $FrR2 $Frd2 $Fy \
						$epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;

nDMaterial NonlinearBS $IDMatSlab $fclU $epslU $fc2U $eps2U $lambda $ftU $Ets $rD1 $SR1 $rdmm1 $rD2 $SR2 $rdmm2 $Fy \
						$epsY $Fy1 $epsY1 $Fu $epsU -$Fy -$epsY -$Fy1 -$epsY1 -$Fu -$epsU $pinchX $pinchY $damage1 $damage2 $betaMUsteel;
#
# Define parameters
set Quad "quad"

if { $Quad == "quad" } {
  set eleArgs "PlaneStress2D $IDMat";
  #set eleArgs "PlaneStrain2D 1"
}


#GEOMETRY
# note : Node number = Element num - 10000

#Element Thickness -Web
set thickness 2.362
set thicknessFlg 2.362
set thicknessSlab 9.84
set SideFlgWidth [expr 23.622*3/8/2*$in]
set WallHeight [expr 47.244*$in]
set WallWidth [expr 23.622*5/8*$in]
set BotSlabHeight [expr 8.86*$in]
#1. Bot slab-Left
set BotSlabLStartX $SideFlgWidth
set BotSlabLStartY 0.
set BotSlabLWidth [expr 23.622*3/8/2*$in]
set BotSlabLHeight $BotSlabHeight
set BotSlabLEleNumX 2.
set BotSlabLEleNumY 3.
set BsLx [expr $BotSlabLWidth/$BotSlabLEleNumX]
set BsLy [expr $BotSlabLHeight/$BotSlabLEleNumY]
#2. Bot slab - Middle
set BotSlabStartX [expr $SideFlgWidth+$BotSlabLWidth]
set BotSlabStartY 0.
set BotSlabWidth $WallWidth
set BotSlabHeight $BotSlabLHeight
set BotSlabEleNumX 5.
set BotSlabEleNumY 3.
set Bsx [expr $BotSlabWidth/$BotSlabEleNumX]
set Bsy [expr $BotSlabHeight/$BotSlabEleNumY]
#3. Bot slab - Right side
set BotSlabRStartX [expr $SideFlgWidth+$BotSlabLWidth+$BotSlabWidth]
set BotSlabRStartY 0.
set BotSlabRWidth $BotSlabLWidth
set BotSlabRHeight $BotSlabHeight
set BotSlabREleNumX 2.
set BotSlabREleNumY 3.
set BsRx [expr $BotSlabRWidth/$BotSlabREleNumX]
set BsRy [expr $BotSlabRHeight/$BotSlabREleNumY]
#1. Top slab - Left side
set TopSlabLStartX $SideFlgWidth
set TopSlabLStartY [expr $BotSlabHeight+$WallHeight]
set TopSlabLWidth [expr 23.622*3/8/2*$in]
set TopSlabLHeight [expr 4.92*$in]
set TopSlabLEleNumX 2.
set TopSlabLEleNumY 2.
set TsLx [expr $TopSlabLWidth/$TopSlabLEleNumX]
set TsLy [expr $TopSlabLHeight/$TopSlabLEleNumY]
#2. Top slab - Middle
set TopSlabStartX [expr $SideFlgWidth+$TopSlabLWidth]
set TopSlabStartY [expr $BotSlabHeight+$WallHeight]
set TopSlabWidth $WallWidth
set TopSlabHeight $TopSlabLHeight
set TopSlabEleNumX 5.
set TopSlabEleNumY 2.
set Tsx [expr $TopSlabWidth/$TopSlabEleNumX]
set Tsy [expr $TopSlabHeight/$TopSlabEleNumY]
#3. Top slab - Right side
set TopSlabRStartX [expr $SideFlgWidth+$TopSlabLWidth+$TopSlabWidth]
set TopSlabRStartY [expr $BotSlabHeight+$WallHeight]
set TopSlabRWidth $TopSlabLWidth
set TopSlabRHeight $TopSlabLHeight
set TopSlabREleNumX 2.
set TopSlabREleNumY 2.
set TsRx [expr $TopSlabRWidth/$TopSlabREleNumX]
set TsRy [expr $TopSlabRHeight/$TopSlabREleNumY]
#4. Flange - Left side
set LeftFlgStartX $SideFlgWidth
set LeftFlgStartY $BotSlabHeight
set LeftFlgWidth $TopSlabLWidth
set LeftFlgHeight $WallHeight
set LeftFlgEleNumX 2.
set LeftFlgEleNumY 12.
set Lfx [expr $LeftFlgWidth/$LeftFlgEleNumX]
set Lfy [expr $LeftFlgHeight/$LeftFlgEleNumY]
#5. Flange - Right side
set RightFlgStartX [expr $SideFlgWidth+$LeftFlgWidth+$WallWidth]
set RightFlgStartY $BotSlabHeight
set RightFlgWidth $TopSlabLWidth
set RightFlgHeight $WallHeight
set RightFlgEleNumX 2.
set RightFlgEleNumY 12.
set Rfx [expr $RightFlgWidth/$RightFlgEleNumX]
set Rfy [expr $RightFlgHeight/$RightFlgEleNumY]
#6. Web or Shear Wall
set WebStartX [expr $SideFlgWidth+$LeftFlgWidth]
set WebStartY $BotSlabHeight
set WebWidth $WallWidth
set WebHeight $WallHeight
set WebEleNumX 5.
set WebEleNumY 12.
set nx [expr $WebWidth/$WebEleNumX]
set ny [expr $WebHeight/$WebEleNumY]
	
set Xmass 2000;
set Ymass 2000;
# NODE DEFINITION
# Top and Bottom Slab
# Top Slab - Sides
# Node number : 2 00(x coord) 00(y coord)
node 201 0 [expr $BotSlabHeight+$WebHeight+0*$Tsy] -mass $Xmass $Ymass
node 202 0 [expr $BotSlabHeight+$WebHeight+1.*$Tsy] -mass $Xmass $Ymass
node 203 0 [expr $BotSlabHeight+$WebHeight+2.*$Tsy] -mass $Xmass $Ymass
node 2101 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr $BotSlabHeight+$WebHeight+0*$Tsy] -mass $Xmass $Ymass
node 2102 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr $BotSlabHeight+$WebHeight+1.*$Tsy] -mass $Xmass $Ymass
node 2103 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr $BotSlabHeight+$WebHeight+2.*$Tsy] -mass $Xmass $Ymass
# Top Slab - Middle
node 80001 [expr $TopSlabLStartX+0.*$TsLx] [expr $TopSlabLStartY+$TsLy] -mass $Xmass $Ymass
node 80101 [expr $TopSlabLStartX+1.*$TsLx] [expr $TopSlabLStartY+$TsLy] -mass $Xmass $Ymass
node 80201 [expr $TopSlabLStartX+2 *$TsLx] [expr $TopSlabLStartY+$TsLy] -mass $Xmass $Ymass
node 80301 [expr $TopSlabStartX+1.*$Tsx] [expr $TopSlabStartY+$Tsy] -mass $Xmass $Ymass
node 80401 [expr $TopSlabStartX+2.*$Tsx] [expr $TopSlabStartY+$Tsy] -mass $Xmass $Ymass
node 80501 [expr $TopSlabStartX+3.*$Tsx] [expr $TopSlabStartY+$Tsy] -mass $Xmass $Ymass
node 80601 [expr $TopSlabStartX+4.*$Tsx] [expr $TopSlabStartY+$Tsy] -mass $Xmass $Ymass
node 80701 [expr $TopSlabRStartX+0.*$TsRx] [expr $TopSlabRStartY+$TsRy] -mass $Xmass $Ymass
node 80801 [expr $TopSlabRStartX+1.*$TsRx] [expr $TopSlabRStartY+$TsRy] -mass $Xmass $Ymass
node 80901 [expr $TopSlabRStartX+2.*$TsRx] [expr $TopSlabRStartY+$TsRy] -mass $Xmass $Ymass
node 80002 [expr $TopSlabLStartX+0.*$TsLx] [expr $TopSlabLStartY+2.*$TsLy] -mass $Xmass $Ymass
node 80102 [expr $TopSlabLStartX+1.*$TsLx] [expr $TopSlabLStartY+2.*$TsLy] -mass $Xmass $Ymass
node 80202 [expr $TopSlabLStartX+2.*$TsLx] [expr $TopSlabLStartY+2.*$TsLy] -mass $Xmass $Ymass
node 80302 [expr $TopSlabStartX+1.*$Tsx] [expr $TopSlabStartY+2.*$Tsy] -mass $Xmass $Ymass
node 80402 [expr $TopSlabStartX+2.*$Tsx] [expr $TopSlabStartY+2.*$Tsy] -mass $Xmass $Ymass
node 80502 [expr $TopSlabStartX+3.*$Tsx] [expr $TopSlabStartY+2.*$Tsy] -mass $Xmass $Ymass
node 80602 [expr $TopSlabStartX+4.*$Tsx] [expr $TopSlabStartY+2.*$Tsy] -mass $Xmass $Ymass
node 80702 [expr $TopSlabRStartX+0.*$TsRx] [expr $TopSlabRStartY+2.*$TsRy] -mass $Xmass $Ymass
node 80802 [expr $TopSlabRStartX+1.*$TsRx] [expr $TopSlabRStartY+2.*$TsRy] -mass $Xmass $Ymass
node 80902 [expr $TopSlabRStartX+2.*$TsRx] [expr $TopSlabRStartY+2.*$TsRy] -mass $Xmass $Ymass
# Bot Slab - Sides
# Node number: 2 00(x coord) 00(y coord)
node 1201 0 [expr 0+0*$Bsy] -mass $Xmass $Ymass
node 1202 0 [expr 0+1.*$Bsy] -mass $Xmass $Ymass
node 1203 0 [expr 0+2.*$Bsy] -mass $Xmass $Ymass
node 1204 0 [expr 0+3.*$Bsy] -mass $Xmass $Ymass
node 12101 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr 0+0*$Bsy] -mass $Xmass $Ymass
node 12102 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr 0+1.*$Bsy] -mass $Xmass $Ymass
node 12103 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr 0+2.*$Bsy] -mass $Xmass $Ymass
node 12104 [expr $SideFlgWidth+$LeftFlgWidth+$WebWidth+$RightFlgWidth+$SideFlgWidth] [expr 0+3.*$Bsy] -mass $Xmass $Ymass
# Bot Slab - Middle
node 100000 [expr $BotSlabLStartX+0.*$BsLx] [expr $BotSlabLStartY] -mass $Xmass $Ymass
node 100100 [expr $BotSlabLStartX+1.*$BsLx] [expr $BotSlabLStartY] -mass $Xmass $Ymass
node 100200 [expr $BotSlabLStartX+2.*$BsLx] [expr $BotSlabLStartY] -mass $Xmass $Ymass
node 100300 [expr $BotSlabStartX+1.*$Bsx] [expr $BotSlabStartY] -mass $Xmass $Ymass
node 100400 [expr $BotSlabStartX+2.*$Bsx] [expr $BotSlabStartY] -mass $Xmass $Ymass
node 100500 [expr $BotSlabStartX+3.*$Bsx] [expr $BotSlabStartY] -mass $Xmass $Ymass
node 100600 [expr $BotSlabStartX+4.*$Bsx] [expr $BotSlabStartY] -mass $Xmass $Ymass
node 100700 [expr $BotSlabRStartX+0.*$BsRx] [expr $BotSlabRStartY] -mass $Xmass $Ymass
node 100800 [expr $BotSlabRStartX+1.*$BsRx] [expr $BotSlabRStartY] -mass $Xmass $Ymass
node 100900 [expr $BotSlabRStartX+2.*$BsRx] [expr $BotSlabRStartY] -mass $Xmass $Ymass
node 100001 [expr $BotSlabLStartX+0.*$BsLx] [expr $BotSlabLStartY+1.*$BsLy] -mass $Xmass $Ymass
node 100101 [expr $BotSlabLStartX+1.*$BsLx] [expr $BotSlabLStartY+1.*$BsLy] -mass $Xmass $Ymass
node 100201 [expr $BotSlabLStartX+2.*$BsLx] [expr $BotSlabLStartY+1.*$BsLy] -mass $Xmass $Ymass
node 100301 [expr $BotSlabStartX+1.*$Bsx] [expr $BotSlabStartY+1.*$Bsy] -mass $Xmass $Ymass
node 100401 [expr $BotSlabStartX+2.*$Bsx] [expr $BotSlabStartY+1.*$Bsy] -mass $Xmass $Ymass
node 100501 [expr $BotSlabStartX+3.*$Bsx] [expr $BotSlabStartY+1.*$Bsy] -mass $Xmass $Ymass
node 100601 [expr $BotSlabStartX+4.*$Bsx] [expr $BotSlabStartY+1.*$Bsy] -mass $Xmass $Ymass
node 100701 [expr $BotSlabRStartX+0.*$BsRx] [expr $BotSlabRStartY+1.*$BsRy] -mass $Xmass $Ymass
node 100801 [expr $BotSlabRStartX+1.*$BsRx] [expr $BotSlabRStartY+1.*$BsRy] -mass $Xmass $Ymass
node 100901 [expr $BotSlabRStartX+2.*$BsRx] [expr $BotSlabRStartY+1.*$BsRy] -mass $Xmass $Ymass
node 100002 [expr $BotSlabLStartX+0.*$BsLx] [expr $BotSlabLStartY+2.*$BsLy] -mass $Xmass $Ymass
node 100102 [expr $BotSlabLStartX+1.*$BsLx] [expr $BotSlabLStartY+2.*$BsLy] -mass $Xmass $Ymass
node 100202 [expr $BotSlabLStartX+2.*$BsLx] [expr $BotSlabLStartY+2.*$BsLy] -mass $Xmass $Ymass
node 100302 [expr $BotSlabStartX+1.*$Bsx] [expr $BotSlabStartY+2.*$Bsy] -mass $Xmass $Ymass
node 100402 [expr $BotSlabStartX+2.*$Bsx] [expr $BotSlabStartY+2.*$Bsy] -mass $Xmass $Ymass
node 100502 [expr $BotSlabStartX+3.*$Bsx] [expr $BotSlabStartY+2.*$Bsy] -mass $Xmass $Ymass
node 100602 [expr $BotSlabStartX+4.*$Bsx] [expr $BotSlabStartY+2.*$Bsy] -mass $Xmass $Ymass
node 100702 [expr $BotSlabRStartX+0.*$BsRx] [expr $BotSlabRStartY+2.*$BsRy] -mass $Xmass $Ymass
node 100802 [expr $BotSlabRStartX+1.*$BsRx] [expr $BotSlabRStartY+2.*$BsRy] -mass $Xmass $Ymass
node 100902 [expr $BotSlabRStartX+2.*$BsRx] [expr $BotSlabRStartY+2.*$BsRy] -mass $Xmass $Ymass

# Left Flange
# Node number : 2 00(x coord) 00(y coord)
node 20000 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+0.]       -mass $Xmass $Ymass
node 20100 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+0.]       -mass $Xmass $Ymass
node 20001 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+$Lfy]     -mass $Xmass $Ymass
node 20101 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+$Lfy]     -mass $Xmass $Ymass
node 20002 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+2.*$Lfy]  -mass $Xmass $Ymass
node 20102 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+2.*$Lfy]  -mass $Xmass $Ymass
node 20003 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+3.*$Lfy]  -mass $Xmass $Ymass
node 20103 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+3.*$Lfy]  -mass $Xmass $Ymass
node 20004 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+4.*$Lfy]  -mass $Xmass $Ymass
node 20104 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+4.*$Lfy]  -mass $Xmass $Ymass
node 20005 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+5.*$Lfy]  -mass $Xmass $Ymass
node 20105 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+5.*$Lfy]  -mass $Xmass $Ymass
node 20006 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+6.*$Lfy]  -mass $Xmass $Ymass
node 20106 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+6.*$Lfy]  -mass $Xmass $Ymass
node 20007 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+7.*$Lfy]  -mass $Xmass $Ymass
node 20107 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+7.*$Lfy]  -mass $Xmass $Ymass
node 20008 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+8.*$Lfy]  -mass $Xmass $Ymass
node 20108 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+8.*$Lfy]  -mass $Xmass $Ymass
node 20009 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+9.*$Lfy]  -mass $Xmass $Ymass
node 20109 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+9.*$Lfy]  -mass $Xmass $Ymass
node 20010 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+10.*$Lfy] -mass $Xmass $Ymass
node 20110 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+10.*$Lfy] -mass $Xmass $Ymass
node 20011 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+11.*$Lfy] -mass $Xmass $Ymass
node 20111 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+11.*$Lfy] -mass $Xmass $Ymass
node 20012 [expr $LeftFlgStartX+0.*$Lfx]  [expr $LeftFlgStartY+12.*$Lfy] -mass $Xmass $Ymass
node 20112 [expr $LeftFlgStartX+1.*$Lfx]  [expr $LeftFlgStartY+12.*$Lfy] -mass $Xmass $Ymass

# Right Flange
# Node number: 6 00(x coord) 00(y coord)
node 60100 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+0.] -mass $Xmass $Ymass
node 60200 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+0.] -mass $Xmass $Ymass
node 60101 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+$Rfy] -mass $Xmass $Ymass
node 60201 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+$Rfy] -mass $Xmass $Ymass
node 60102 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+2.*$Rfy] -mass $Xmass $Ymass
node 60202 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+2.*$Rfy] -mass $Xmass $Ymass
node 60103 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+3.*$Rfy] -mass $Xmass $Ymass
node 60203 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+3.*$Rfy] -mass $Xmass $Ymass
node 60104 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+4.*$Rfy] -mass $Xmass $Ymass
node 60204 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+4.*$Rfy] -mass $Xmass $Ymass
node 60105 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+5.*$Rfy] -mass $Xmass $Ymass
node 60205 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+5.*$Rfy] -mass $Xmass $Ymass
node 60106 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+6.*$Rfy] -mass $Xmass $Ymass
node 60206 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+6.*$Rfy] -mass $Xmass $Ymass
node 60107 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+7.*$Rfy] -mass $Xmass $Ymass
node 60207 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+7.*$Rfy] -mass $Xmass $Ymass
node 60108 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+8.*$Rfy] -mass $Xmass $Ymass
node 60208 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+8.*$Rfy] -mass $Xmass $Ymass
node 60109 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+9.*$Rfy] -mass $Xmass $Ymass
node 60209 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+9.*$Rfy] -mass $Xmass $Ymass
node 60110 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+10.*$Rfy] -mass $Xmass $Ymass
node 60210 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+10.*$Rfy] -mass $Xmass $Ymass
node 60111 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+11.*$Rfy] -mass $Xmass $Ymass
node 60211 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+11.*$Rfy] -mass $Xmass $Ymass
node 60112 [expr $RightFlgStartX+1.*$Rfx] [expr $RightFlgStartY+12.*$Rfy] -mass $Xmass $Ymass
node 60212 [expr $RightFlgStartX+2.*$Rfx] [expr $RightFlgStartY+12.*$Rfy] -mass $Xmass $Ymass
# Web - Node number : 4 00(x coord) 00(y coord)
node 40000 [expr $WebStartX+0.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40100 [expr $WebStartX+1.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40200 [expr $WebStartX+2.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40300 [expr $WebStartX+3.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40400 [expr $WebStartX+4.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40500 [expr $WebStartX+5.*$nx] [expr $WebStartY+0] -mass $Xmass $Ymass
node 40001 [expr $WebStartX+0.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40101 [expr $WebStartX+1.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40201 [expr $WebStartX+2.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40301 [expr $WebStartX+3.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40401 [expr $WebStartX+4.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40501 [expr $WebStartX+5.*$nx] [expr $WebStartY+$ny] -mass $Xmass $Ymass
node 40002 [expr $WebStartX+0.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40102 [expr $WebStartX+1.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40202 [expr $WebStartX+2.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40302 [expr $WebStartX+3.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40402 [expr $WebStartX+4.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40502 [expr $WebStartX+5.*$nx] [expr $WebStartY+2.*$ny] -mass $Xmass $Ymass
node 40003 [expr $WebStartX+0.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40103 [expr $WebStartX+1.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40203 [expr $WebStartX+2.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40303 [expr $WebStartX+3.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40403 [expr $WebStartX+4.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40503 [expr $WebStartX+5.*$nx] [expr $WebStartY+3.*$ny] -mass $Xmass $Ymass
node 40004 [expr $WebStartX+0.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40104 [expr $WebStartX+1.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40204 [expr $WebStartX+2.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40304 [expr $WebStartX+3.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40404 [expr $WebStartX+4.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40504 [expr $WebStartX+5.*$nx] [expr $WebStartY+4.*$ny] -mass $Xmass $Ymass
node 40005 [expr $WebStartX+0.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40105 [expr $WebStartX+1.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40205 [expr $WebStartX+2.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40305 [expr $WebStartX+3.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40405 [expr $WebStartX+4.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40505 [expr $WebStartX+5.*$nx] [expr $WebStartY+5.*$ny] -mass $Xmass $Ymass
node 40006 [expr $WebStartX+0.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40106 [expr $WebStartX+1.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40206 [expr $WebStartX+2.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40306 [expr $WebStartX+3.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40406 [expr $WebStartX+4.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40506 [expr $WebStartX+5.*$nx] [expr $WebStartY+6.*$ny] -mass $Xmass $Ymass
node 40007 [expr $WebStartX+0.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40107 [expr $WebStartX+1.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40207 [expr $WebStartX+2.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40307 [expr $WebStartX+3.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40407 [expr $WebStartX+4.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40507 [expr $WebStartX+5.*$nx] [expr $WebStartY+7.*$ny] -mass $Xmass $Ymass
node 40008 [expr $WebStartX+0.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40108 [expr $WebStartX+1.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40208 [expr $WebStartX+2.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40308 [expr $WebStartX+3.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40408 [expr $WebStartX+4.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40508 [expr $WebStartX+5.*$nx] [expr $WebStartY+8.*$ny] -mass $Xmass $Ymass
node 40009 [expr $WebStartX+0.*$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40109 [expr $WebStartX+1.*$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40209 [expr $WebStartX+2.*$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40309 [expr $WebStartX+3 *$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40409 [expr $WebStartX+4.*$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40509 [expr $WebStartX+5.*$nx] [expr $WebStartY+9.*$ny] -mass $Xmass $Ymass
node 40010 [expr $WebStartX+0.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40110 [expr $WebStartX+1.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40210 [expr $WebStartX+2.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40310 [expr $WebStartX+3.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40410 [expr $WebStartX+4.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40510 [expr $WebStartX+5.*$nx] [expr $WebStartY+10.*$ny] -mass $Xmass $Ymass
node 40011 [expr $WebStartX+0.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40111 [expr $WebStartX+1.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40211 [expr $WebStartX+2.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40311 [expr $WebStartX+3.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40411 [expr $WebStartX+4.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40511 [expr $WebStartX+5.*$nx] [expr $WebStartY+11.*$ny] -mass $Xmass $Ymass
node 40012 [expr $WebStartX+0.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
node 40112 [expr $WebStartX+1.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
node 40212 [expr $WebStartX+2.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
node 40312 [expr $WebStartX+3.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
node 40412 [expr $WebStartX+4.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
node 40512 [expr $WebStartX+5.*$nx] [expr $WebStartY+12.*$ny] -mass $Xmass $Ymass
# ELEMENT DEFINITION
# <Bot slab>
# Bot slab-Side
#element quad SeleTag SiNode SjNode SkNode SiNode Sthick Stype SmatTag
element quad 1301 1201 100000 100001 1202 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 1302 1202 100001 100002 1203 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 1303 1203 100002 20000 1204 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 13101 100900 12101 12102 100901 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 13102 100901 12102 12103 100902 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 13103 100902 12103 12104 60200 $thicknessSlab "PlaneStress2D" $IDMatSlab
# Bot slab - Middle
element quad 110101 100000 100100 100101 100001 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110201 100100 100200 100201 100101 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110301 100200 100300 100301 100201 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110401 100300 100400 100401 100301 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110501 100400 100500 100501 100401 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110601 100500 100600 100601 100501 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110701 100600 100700 100701 100601 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110801 100700 100800 100801 100701 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110901 100800 100900 100901 100801 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110102 100001 100101 100102 100002 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110202 100101 100201 100202 100102 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110302 100201 100301 100302 100202 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110402 100301 100401 100402 100302 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110502 100401 100501 100502 100402 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110602 100501 100601 100602 100502 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110702 100601 100701 100702 100602 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110802 100701 100801 100802 100702 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110902 100801 100901 100902 100802 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 110103 100002 100102 20100 20000 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110203 100102 100202 40000 20100 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110303 100202 100302 40100 40000 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110403 100302 100402 40200 40100 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110503 100402 100502 40300 40200 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110603 100502 100602 40400 40300 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110703 100602 100702 40500 40400 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110803 100702 100802 60100 40500 $thicknessSlab "PlaneStress2D" $IDMatSlab  
element quad 110903 100802 100902 60200 60100 $thicknessSlab "PlaneStress2D" $IDMatSlab  
# <Top slab>
#element quad SeleTag SiNode SjNode SkNode SlNode Sthick Stype SmatTag
#-Side Slab
element quad 301 201 20012 80001 202 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 302 202 80001 80002 203 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 3101 60212 2101 2102 80901 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 3102 80901 2102 2103 80902 $thicknessSlab "PlaneStress2D" $IDMatSlab
#-Middle Slab
element quad 90101 20012 20112 80101 80001 $thicknessSlab "PlaneStress2D" $IDMatSlab
element quad 90201 20112 40012 80201 80101 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90301 40012 40112 80301 80201 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90401 40112 40212 80401 80301 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90501 40212 40312 80501 80401 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90601 40312 40412 80601 80501 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90701 40412 40512 80701 80601 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90801 40512 60112 80801 80701 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90901 60112 60212 80901 80801 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90102 80001 80101 80102 80002 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90202 80101 80201 80202 80102 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90302 80201 80301 80302 80202 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90402 80301 80401 80402 80302 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90502 80401 80501 80502 80402 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90602 80501 80601 80602 80502 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90702 80601 80701 80702 80602 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90802 80701 80801 80802 80702 $thicknessSlab "PlaneStress2D" $IDMatSlab 
element quad 90902 80801 80901 80902 80802 $thicknessSlab "PlaneStress2D" $IDMatSlab 
# <Left Flange>
element quad 30101 20000 20100 20101 20001 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30201 20100 40000 40001 20101 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30102 20001 20101 20102 20002 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30202 20101 40001 40002 20102 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30103 20002 20102 20103 20003 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30203 20102 40002 40003 20103 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30104 20003 20103 20104 20004 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30204 20103 40003 40004 20104 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30105 20004 20104 20105 20005 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30205 20104 40004 40005 20105 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30106 20005 20105 20106 20006 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30206 20105 40005 40006 20106 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30107 20006 20106 20107 20007 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30207 20106 40006 40007 20107 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30108 20007 20107 20108 20008 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30208 20107 40007 40008 20108 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30109 20008 20108 20109 20009 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30209 20108 40008 40009 20109 $thicknessFlg "PlaneStress2D" $IDMatFlg   
element quad 30110 20009 20109 20110 20010 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 30210 20109 40009 40010 20110 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 30111 20010 20110 20111 20011 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 30211 20110 40010 40011 20111 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 30112 20011 20111 20112 20012 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 30212 20111 40011 40012 20112 $thicknessFlg "PlaneStress2D" $IDMatFlg
# <Right Flange>
element quad 70101 40500 60100 60101 40501 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70201 60100 60200 60201 60101 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70102 40501 60101 60102 40502 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70202 60101 60201 60202 60102 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70103 40502 60102 60103 40503 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70203 60102 60202 60203 60103 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70104 40503 60103 60104 40504 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70204 60103 60203 60204 60104 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70105 40504 60104 60105 40505 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70205 60104 60204 60205 60105 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70106 40505 60105 60106 40506 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70206 60105 60205 60206 60106 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70107 40506 60106 60107 40507 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70207 60106 60206 60207 60107 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70108 40507 60107 60108 40508 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70208 60107 60207 60208 60108 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70109 40508 60108 60109 40509 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70209 60108 60208 60209 60109 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70110 40509 60109 60110 40510 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70210 60109 60209 60210 60110 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70111 40510 60110 60111 40511 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70211 60110 60210 60211 60111 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70112 40511 60111 60112 40512 $thicknessFlg "PlaneStress2D" $IDMatFlg
element quad 70212 60111 60211 60212 60112 $thicknessFlg "PlaneStress2D" $IDMatFlg
#<Web>
#element quad SeleTag SiNode SjNode SkNode SlNode Sthick Stype SmatTag
element quad 50101 40000 40100 40101 40001 $thickness "PlaneStress2D" $IDMat
element quad 50201 40100 40200 40201 40101 $thickness "PlaneStress2D" $IDMat
element quad 50301 40200 40300 40301 40201 $thickness "PlaneStress2D" $IDMat
element quad 50401 40300 40400 40401 40301 $thickness "PlaneStress2D" $IDMat
element quad 50501 40400 40500 40501 40401 $thickness "PlaneStress2D" $IDMat
element quad 50102 40001 40101 40102 40002 $thickness "PlaneStress2D" $IDMat
element quad 50202 40101 40201 40202 40102 $thickness "PlaneStress2D" $IDMat
element quad 50302 40201 40301 40302 40202 $thickness "PlaneStress2D" $IDMat
element quad 50402 40301 40401 40402 40302 $thickness "PlaneStress2D" $IDMat
element quad 50502 40401 40501 40502 40402 $thickness "PlaneStress2D" $IDMat
element quad 50103 40002 40102 40103 40003 $thickness "PlaneStress2D" $IDMat
element quad 50203 40102 40202 40203 40103 $thickness "PlaneStress2D" $IDMat
element quad 50303 40202 40302 40303 40203 $thickness "PlaneStress2D" $IDMat
element quad 50403 40302 40402 40403 40303 $thickness "PlaneStress2D" $IDMat
element quad 50503 40402 40502 40503 40403 $thickness "PlaneStress2D" $IDMat
element quad 50104 40003 40103 40104 40004 $thickness "PlaneStress2D" $IDMat
element quad 50204 40103 40203 40204 40104 $thickness "PlaneStress2D" $IDMat
element quad 50304 40203 40303 40304 40204 $thickness "PlaneStress2D" $IDMat
element quad 50404 40303 40403 40404 40304 $thickness "PlaneStress2D" $IDMat
element quad 50504 40403 40503 40504 40404 $thickness "PlaneStress2D" $IDMat
element quad 50105 40004 40104 40105 40005 $thickness "PlaneStress2D" $IDMat
element quad 50205 40104 40204 40205 40105 $thickness "PlaneStress2D" $IDMat
element quad 50305 40204 40304 40305 40205 $thickness "PlaneStress2D" $IDMat
element quad 50405 40304 40404 40405 40305 $thickness "PlaneStress2D" $IDMat
element quad 50505 40404 40504 40505 40405 $thickness "PlaneStress2D" $IDMat
element quad 50106 40005 40105 40106 40006 $thickness "PlaneStress2D" $IDMat
element quad 50206 40105 40205 40206 40106 $thickness "PlaneStress2D" $IDMat
element quad 50306 40205 40305 40306 40206 $thickness "PlaneStress2D" $IDMat
element quad 50406 40305 40405 40406 40306 $thickness "PlaneStress2D" $IDMat
element quad 50506 40405 40505 40506 40406 $thickness "PlaneStress2D" $IDMat
element quad 50107 40006 40106 40107 40007 $thickness "PlaneStress2D" $IDMat
element quad 50207 40106 40206 40207 40107 $thickness "PlaneStress2D" $IDMat
element quad 50307 40206 40306 40307 40207 $thickness "PlaneStress2D" $IDMat
element quad 50407 40306 40406 40407 40307 $thickness "PlaneStress2D" $IDMat
element quad 50507 40406 40506 40507 40407 $thickness "PlaneStress2D" $IDMat
element quad 50108 40007 40107 40108 40008 $thickness "PlaneStress2D" $IDMat
element quad 50208 40107 40207 40208 40108 $thickness "PlaneStress2D" $IDMat
element quad 50308 40207 40307 40308 40208 $thickness "PlaneStress2D" $IDMat
element quad 50408 40307 40407 40408 40308 $thickness "PlaneStress2D" $IDMat
element quad 50508 40407 40507 40508 40408 $thickness "PlaneStress2D" $IDMat
element quad 50109 40008 40108 40109 40009 $thickness "PlaneStress2D" $IDMat
element quad 50209 40108 40208 40209 40109 $thickness "PlaneStress2D" $IDMat
element quad 50309 40208 40308 40309 40209 $thickness "PlaneStress2D" $IDMat
element quad 50409 40308 40408 40409 40309 $thickness "PlaneStress2D" $IDMat
element quad 50509 40408 40508 40509 40409 $thickness "PlaneStress2D" $IDMat
element quad 50110 40009 40109 40110 40010 $thickness "PlaneStress2D" $IDMat
element quad 50210 40109 40209 40210 40110 $thickness "PlaneStress2D" $IDMat
element quad 50310 40209 40309 40310 40210 $thickness "PlaneStress2D" $IDMat
element quad 50410 40309 40409 40410 40310 $thickness "PlaneStress2D" $IDMat
element quad 50510 40409 40509 40510 40410 $thickness "PlaneStress2D" $IDMat
element quad 50111 40010 40110 40111 40011 $thickness "PlaneStress2D" $IDMat
element quad 50211 40110 40210 40211 40111 $thickness "PlaneStress2D" $IDMat
element quad 50311 40210 40310 40311 40211 $thickness "PlaneStress2D" $IDMat
element quad 50411 40310 40410 40411 40311 $thickness "PlaneStress2D" $IDMat
element quad 50511 40410 40510 40511 40411 $thickness "PlaneStress2D" $IDMat
element quad 50112 40011 40111 40112 40012 $thickness "PlaneStress2D" $IDMat
element quad 50212 40111 40211 40212 40112 $thickness "PlaneStress2D" $IDMat
element quad 50312 40211 40311 40312 40212 $thickness "PlaneStress2D" $IDMat
element quad 50412 40311 40411 40412 40312 $thickness "PlaneStress2D" $IDMat
element quad 50512 40411 40511 40512 40412 $thickness "PlaneStress2D" $IDMat
# BOUNDARY CONDITION (Single point constraints)
#  node ul u2
fix 100000 1 1
fix 100100 1 1
fix 100200 1 1
fix 100300 1 1
fix 100400 1 1
fix 100500 1 1
fix 100600 1 1
fix 100700 1 1
fix 100800 1 1
fix 100900 1 1
fix 1201 1 1
fix 12101 1 1

puts "Model Built"

## create the display
set displayType "PERSPECTIVE"
recorder display g3 10 10 800 600 -wipe
if {$displayType == "PERSPECTIVE"} {
  prp -12500 -25000 30000
  #vrp 0 -500 250
  vup 0 0 1
  #vpn -1 -1 0.5
  viewWindow -2000 2000 -2000 2000
}
port -1 1 -1 1
projection 1
fill 0
display 1 1 10

#
# Apply the displacement history and Perform a transient analysis
#
set timeIncr 0.00004;
set displacement1 "Series -dt $timeIncr -filePath $DHfile1 -factor 100.01";
pattern MultipleSupport 2 {
	groundMotion 2 Plain -disp $displacement1
	# « Middle eight nodes of the top slab - horizontal »
	#$nodeTag Sdirn SgMotionTag
	imposedMotion 80102 1 2
	imposedMotion 80202 1 2
	imposedMotion 80302 1 2
	imposedMotion 80402 1 2
	imposedMotion 80502 1 2
	imposedMotion 80602 1 2
	imposedMotion 80702 1 2
	imposedMotion 80802 1 2
}

# Analyze
rayleigh 0. 0. 0. 0.05
constraints Penalty 1.0e12 1.0e12
numberer RCM
system UmfPack
# Convergence test
# tolerance maxlter displayCode
test EnergyIncr 1.0e-10 20 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient
#Pushover analysis or
#analyze 1901 Stimelncr
#Cyclic analysis
analyze 10161 $timeIncr