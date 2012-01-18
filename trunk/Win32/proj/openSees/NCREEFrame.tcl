# NCREE 3D Frame
# applied axial load first

# ------------------------------
# Units: N, mm, sec, MPa
# ------------------------------

# ------------------------------
# Start of model generation
# ------------------------------



# ------------------------------------------------------------------------------
# Create ModelBuilder (with three-dimensions and 6 dof/node)
# ------------------------------------------------------------------------------
  model basic -ndm 3 -ndf 6

# ---------------------------------------------------------
# Create nodes (6dof)
# ---------------------------------------------------------

# Set the dimensions of the frame

  set deltaAB 3500;   # Distance between Column Lines A and B
  set deltaBC 5500;   # Distance between Column Lines B and C
  set delta12 3000;   # Distance between Column Lines 1 and 2
  set cH1 2000;      # Height to bottom of 2nd floor slab
  set cH2 2200;      # Height to center of 2nd floor slab
  set cH3 2400;      # Height to top of 2nd floor slab
  set cH4 3900;      # Height to bottom of 3rd floor slab
  set cH5 4100;      # Height to center of 3rd floor slab
  set cH6 4700;      # Height to top of columns

# Set locations for load placement

  set x1 750;
  set x2 1750;
  set x3 3250;
  set x4 7600;
  set x5 8500;
  set y1 600;
  set y2 1500;
  set y3 2400;

# Create the nodes (one node at the center of the end of each member)
# Nodes 101 < node < 166 represent connections used for nonlinear beamcolumn members
# Nodes 167 < node < 198 represent nodes used for loads

  node 101 0 0 0;
  node 102 $deltaAB 0 0;
  node 103 [expr $deltaAB+$deltaBC] 0 0;
  node 104 0 $delta12 0;
  node 105 $deltaAB $delta12 0;
  node 106 [expr $deltaAB+$deltaBC] $delta12 0;
  node 111 0 0 $cH1;
  node 112 $deltaAB 0 $cH1;
  node 113 [expr $deltaAB+$deltaBC] 0 $cH1;
  node 114 0 $delta12 $cH1;
  node 115 $deltaAB $delta12 $cH1;
  node 116 [expr $deltaAB+$deltaBC] $delta12 $cH1;
  node 121 0 0 $cH2;
  node 122 $deltaAB 0 $cH2;
  node 123 [expr $deltaAB+$deltaBC] 0 $cH2;
  node 124 0 $delta12 $cH2;
  node 125 $deltaAB $delta12 $cH2;
  node 126 [expr $deltaAB+$deltaBC] $delta12 $cH2;
  node 131 0 0 $cH3;
  node 132 $deltaAB 0 $cH3;
  node 133 [expr $deltaAB+$deltaBC] 0 $cH3;
  node 134 0 $delta12 $cH3;
  node 135 $deltaAB $delta12 $cH3;
  node 136 [expr $deltaAB+$deltaBC] $delta12 $cH3;
  node 141 0 0 $cH4;
  node 142 $deltaAB 0 $cH4;
  node 143 [expr $deltaAB+$deltaBC] 0 $cH4;
  node 144 0 $delta12 $cH4;
  node 145 $deltaAB $delta12 $cH4;
  node 146 [expr $deltaAB+$deltaBC] $delta12 $cH4;
  node 151 0 0 $cH5;
  node 152 $deltaAB 0 $cH5;
  node 153 [expr $deltaAB+$deltaBC] 0 $cH5;
  node 154 0  $delta12 $cH5;
  node 155 $deltaAB $delta12 $cH5 $delta12;
  node 156 [expr $deltaAB+$deltaBC] $delta12 $cH5;
  node 161 0 0 $cH6;
  node 162 $deltaAB 0 $cH6;
  node 163 [expr $deltaAB+$deltaBC] 0 $cH6;
  node 164 0  $delta12 $cH6;
  node 165 $deltaAB $delta12 $cH6;
  node 166 [expr $deltaAB+$deltaBC] $delta12 $cH6;

  node 167 $x1 0 $cH2;
  node 168 $x2 0 $cH2;
  node 169 $x3 0 $cH2;
  node 170 $x4 0 $cH2;
  node 171 $x5 0 $cH2;
  node 172 $x1 $delta12 $cH2;
  node 173 $x2 $delta12 $cH2;
  node 174 $x3 $delta12 $cH2;
  node 175 $x4 $delta12 $cH2;
  node 176 $x5 $delta12 $cH2;
  node 177 0 $y1 $cH2;
  node 178 0 $y2 $cH2;
  node 179 0 $y3 $cH2;
  node 180 $deltaAB $y1 $cH2;
  node 181 $deltaAB $y2 $cH2;
  node 182 $deltaAB $y3 $cH2;
  node 183 $x1 0 $cH5;
  node 184 $x2 0 $cH5;
  node 185 $x3 0 $cH5;
  node 186 $x4 0 $cH5;
  node 187 $x5 0 $cH5;
  node 188 $x1 $delta12 $cH5;
  node 189 $x2 $delta12 $cH5;
  node 190 $x3 $delta12 $cH5;
  node 191 $x4 $delta12 $cH5;
  node 192 $x5 $delta12 $cH5;
  node 193 0 $y1 $cH5;
  node 194 0 $y2 $cH5;
  node 195 0 $y3 $cH5;
  node 196 $deltaAB $y1 $cH5;
  node 197 $deltaAB $y2 $cH5;
  node 198 $deltaAB $y3 $cH5;



# Fix supports and connections

  fix 101 1 1 1 1 1 1;    # fully-fixed
  fix 102 1 1 1 1 1 1;
  fix 103 1 1 1 1 1 1;
  fix 104 1 1 1 1 1 1;
  fix 105 1 1 1 1 1 1;
  fix 106 1 1 1 1 1 1;
  

# --------------------------------------------------
# Define Materials for Nonlinear Members
# --------------------------------------------------

# CONCRETE tag f'c ec0 f'cu ecu

# 1 for unconfined concrete fc'= 27.5
# 2 for confined concrete, fc' = 27.5*variable dependant on tangental reinforcement

  set fc 27.5;       # unconfined concrete
  set fc80 [expr 1.20*$fc];   # confined concrete in columns, #3 @ 80  !!!VARIABLE NEEDS TO BE DETERMINED!!!
  set fc100 [expr 1.15*$fc];   # confined concrete in columns, #3 @ 100  !!!VARIABLE NEEDS TO BE DETERMINED!!!
  set fc280 [expr 1.10*$fc];   # confined concrete in columns, #3 @ 100  !!!VARIABLE NEEDS TO BE DETERMINED!!!
  set fcw [expr 1.10*$fc];  # confined concrete in WallS  !!!VARIABLE NEEDS TO BE DETERMINED!!!


# Cover concrete (unconfined)
  uniaxialMaterial Concrete01 1 [expr -$fc] -0.003  -6.4 -0.006

# Core concrete (confined #3@80) 
  uniaxialMaterial Concrete01 2 [expr -$fc80] -0.003 -7.0 -0.006

# Core concrete (confined #3@100) 
  uniaxialMaterial Concrete01 3 [expr -$fc100] -0.003 -7.0 -0.006

# Core concrete (confined #3@280) 
  uniaxialMaterial Concrete01 4 [expr -$fc280] -0.003 -7.0 -0.006

# STEEL

  set fy 411.9; # Yield stress 
  set E 200000.0; # Young's modulus

# tag fy E0 b
#  uniaxialMaterial Steel01 3 $fy $E 0.024

#                               tag  fy   E0  fpc  rou
     uniaxialMaterial  SteelZ01   5   $fy  $E  32.0 0.017

#   uniaxialMaterial  Steel01    5   [expr $fy]  $E  0.02


# -----------------------------------------------
# Define Fibers
# -----------------------------------------------

# Columns
# -----------------------------------------------

# Column paramaters

  set cWidth 400.0;   # column width
  set cDepth 400.0;   # column height
  set cover 50.0;   # cover
  set As6 285.0;    # area of no. 6 bars

# Variables derived from the parameters

  set cy1 [expr $cDepth/2.0];
  set cz1 [expr $cWidth/2.0];

  # Create a Uniaxial Fiber object
  section Fiber 1 {   # #3@80

  # Create the concrete core fibers
  #         mat num num
  patch rect 2  2  2 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover];

  # Create the concrete cover fibers (top, bottom, left, right)
  patch rect 1 2 1 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1;
  patch rect 1 2 1 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1];
  patch rect 1 2 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover];
  patch rect 1 2 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover];

  # Create the reinforcing fibers (4 layers)
  layer straight 5 4 $As6 [expr $cy1-$cover] [expr $cz1-$cover] [expr $cy1-$cover] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr -($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr -($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr ($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr ($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 4 $As6 [expr $cover-$cy1] [expr $cz1-$cover] [expr $cover-$cy1] [expr $cover-$cz1];

  } 

  # Create a Uniaxial Fiber object
  section Fiber 2 {   # #3@100

  # Create the concrete core fibers
  #         mat num num
  patch rect 3  2  2 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover];

  # Create the concrete cover fibers (top, bottom, left, right)
  patch rect 1 2 1 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1;
  patch rect 1 2 1 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1];
  patch rect 1 2 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover];
  patch rect 1 2 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover];

  # Create the reinforcing fibers (4 layers)
  layer straight 5 4 $As6 [expr $cy1-$cover] [expr $cz1-$cover] [expr $cy1-$cover] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr -($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr -($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr ($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr ($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 4 $As6 [expr $cover-$cy1] [expr $cz1-$cover] [expr $cover-$cy1] [expr $cover-$cz1];

  } 

  # Create a Uniaxial Fiber object
  section Fiber 3 {   # #3@280

  # Create the concrete core fibers
  #         mat num num
  patch rect 4  2  2 [expr $cover-$cy1] [expr $cover-$cz1] [expr $cy1-$cover] [expr $cz1-$cover];

  # Create the concrete cover fibers (top, bottom, left, right)
  patch rect 1 2 1 [expr -$cy1] [expr $cz1-$cover] $cy1 $cz1;
  patch rect 1 2 1 [expr -$cy1] [expr -$cz1] $cy1 [expr $cover-$cz1];
  patch rect 1 2 1 [expr -$cy1] [expr $cover-$cz1] [expr $cover-$cy1] [expr $cz1-$cover];
  patch rect 1 2 1 [expr $cy1-$cover] [expr $cover-$cz1] $cy1 [expr $cz1-$cover];

  # Create the reinforcing fibers (4 layers)
  layer straight 5 4 $As6 [expr $cy1-$cover] [expr $cz1-$cover] [expr $cy1-$cover] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr -($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr -($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 2 $As6 [expr ($cWidth-2*$cover)/6] [expr $cz1-$cover] [expr ($cWidth-2*$cover)/6] [expr $cover-$cz1];
  layer straight 5 4 $As6 [expr $cover-$cy1] [expr $cz1-$cover] [expr $cover-$cy1] [expr $cover-$cz1];

  } 

# Beams
# -----------------------------------------------

# Beam paramaters

  set bWidth 300.0;   # beam width
  set bDepth 400.0;   # beam height


# Variables derived from the parameters

  set by1 [expr $bDepth/2.0];
  set bz1 [expr $bWidth/2.0];

  # Create a Uniaxial Fiber object
  section Fiber 4 {

  # Create the concrete core fibers
  #         mat num num
  patch rect 3  2  2 [expr $cover-$by1] [expr $cover-$bz1] [expr $by1-$cover] [expr $bz1-$cover];

  # Create the concrete cover fibers (top, bottom, left, right)
  patch rect 1 2 1 [expr -$by1] [expr $bz1-$cover] $by1 $bz1;
  patch rect 1 2 1 [expr -$by1] [expr -$bz1] $by1 [expr $cover-$bz1];
  patch rect 1 2 1 [expr -$by1] [expr $cover-$bz1] [expr $cover-$by1] [expr $bz1-$cover];
  patch rect 1 2 1 [expr $by1-$cover] [expr $cover-$bz1] $by1 [expr $bz1-$cover];

  # Create the reinforcing fibers (4 layers)
  layer straight 3 2 $As6 [expr $by1-$cover] [expr $bz1-$cover] [expr $by1-$cover] [expr $cover-$bz1];
  layer straight 3 2 $As6 [expr -($bWidth-2*$cover)/6] [expr $bz1-$cover] [expr -($bWidth-2*$cover)/6] [expr $cover-$bz1];
  layer straight 3 2 $As6 [expr ($bWidth-2*$cover)/6] [expr $bz1-$cover] [expr ($bWidth-2*$cover)/6] [expr $cover-$bz1];
  layer straight 3 2 $As6 [expr $cover-$by1] [expr $bz1-$cover] [expr $cover-$by1] [expr $cover-$bz1];

  } 

# ----------------------------------------------------------
# Define Section Aggregators
# ---------------------------------------------------------

  set U 1.e10;       # Really large number
  set G $U;      # Torsional stiffness modulus
  set J 1.;      # Torsional stiffness of section
  set GJ [expr $G*$J];   # Torsional stiffness 

  uniaxialMaterial Elastic 6 $GJ;

  section Aggregator 5 6 T -section 1;   # #3@80
  section Aggregator 6 6 T -section 2;    # #3@100
  section Aggregator 7 6 T -section 3;  # #3@280
  section Aggregator 8 6 T -section 4;  # Beam    

# ----------------------------------------------------------
# Define Column Elements
# ---------------------------------------------------------

# Geometry of column elements
   
#                  tag vecxzX vecxzY vecxzZ
  geomTransf Linear 1    0      1       0;

# Number of integration points along length of element

  set np 10;
  set iterNum 10;
  set iterTol 1.e-3;

# Create the beam elements using beam-column elements
# Elements 0 < element < 100 are used for columns
# tag ndI ndJ nsecs secID transfTag

  set j 1;
  while {$j < 7} {
    element nonlinearBeamColumn   $j           [expr 100+$j] [expr 110+$j] $np  6  1  -iter $iterNum $iterTol;
    element nonlinearBeamColumn   [expr 10+$j] [expr 110+$j] [expr 120+$j] $np  5  1  -iter $iterNum $iterTol;
    element nonlinearBeamColumn   [expr 20+$j] [expr 120+$j] [expr 130+$j] $np  5  1  -iter $iterNum $iterTol;
    if {$j == 3} {
      element nonlinearBeamColumn   [expr 30+$j] [expr 130+$j] [expr 140+$j] $np  7  1  -iter $iterNum $iterTol;
    } elseif {$j == 6} {
        element nonlinearBeamColumn   [expr 30+$j] [expr 130+$j] [expr 140+$j] $np  7  1  -iter $iterNum $iterTol;
    } else {
        element nonlinearBeamColumn   [expr 30+$j] [expr 130+$j] [expr 140+$j] $np  6  1  -iter $iterNum $iterTol;
    }
    element nonlinearBeamColumn   [expr 40+$j] [expr 140+$j] [expr 150+$j] $np  5  1  -iter $iterNum $iterTol;
    element nonlinearBeamColumn   [expr 50+$j] [expr 150+$j] [expr 160+$j] $np  5  1  -iter $iterNum $iterTol;
  set j [expr $j +1]
  }


# ----------------------------------------------------------
# Define Beam Elements
# ---------------------------------------------------------

# Geometry of beam elements
# tag 

  geomTransf Linear 2 0 0 -1;


# Create the beam elements using beam-column elements
# Elements 100 < element < 300 are used for beams
# tag ndI ndJ nsecs secID transfTag

  element nonlinearBeamColumn 121 121 167 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 167 167 168 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 168 168 169 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 169 169 122 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 122 122 170 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 170 170 171 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 171 171 123 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 124 124 172 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 172 172 173 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 173 173 174 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 174 174 125 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 125 125 175 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 175 175 176 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 176 176 126 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 151 151 183 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 183 183 184 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 184 184 185 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 185 185 152 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 152 152 186 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 186 186 187 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 187 187 153 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 154 154 188 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 188 188 189 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 189 189 190 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 190 190 155 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 155 155 191 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 191 191 192 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 192 192 156 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 221 121 177 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 277 177 178 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 278 178 179 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 279 179 124 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 222 122 180 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 280 180 181 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 281 181 182 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 282 182 125 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 223 123 126 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 251 151 193 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 293 193 194 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 294 194 195 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 295 195 154 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 252 152 196 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 296 196 197 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 297 197 198 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 298 198 155 $np 8 2 -iter $iterNum $iterTol;
  element nonlinearBeamColumn 253 153 156 $np 8 2 -iter $iterNum $iterTol;


# ---------------------------
# Define horizontal loads
# ---------------------------


  set P3 1000.0;        # Load on 3rd Story
  set P2 [expr $P3/1.86];   # Load on 2nd story

# Create a Plain load pattern with a linear TimeSeries
  pattern Plain 1 "Linear" {

    # Create the nodal load - command: load nodeID xForce yForce zForce
    load 167 0 [expr  5*$P2/48] 0 0 0 0;
    load 168 0 [expr  5*$P2/48] 0 0 0 0; 
    load 169 0 [expr  5*$P2/48] 0 0 0 0;  
    load 170 0 [expr  3*$P2/32] 0 0 0 0;  
    load 171 0 [expr  3*$P2/32] 0 0 0 0;   
    load 172 0 [expr  5*$P2/48] 0 0 0 0; 
    load 173 0 [expr  5*$P2/48] 0 0 0 0;  
    load 174 0 [expr  5*$P2/48] 0 0 0 0;   
    load 175 0 [expr  3*$P2/32] 0 0 0 0;  
    load 176 0 [expr  3*$P2/32] 0 0 0 0; 
    load 177 [expr $P2/6] 0 0 0 0 0;   
    load 178 [expr $P2/6] 0 0 0 0 0; 
    load 179 [expr $P2/6] 0 0 0 0 0;   
    load 180 [expr $P2/6] 0 0 0 0 0; 
    load 181 [expr $P2/6] 0 0 0 0 0;   
    load 182 [expr $P2/6] 0 0 0 0 0;   
    load 183 0 [expr  5*$P3/48] 0 0 0 0; 
    load 184 0 [expr  5*$P3/48] 0 0 0 0;  
    load 185 0 [expr  5*$P3/48] 0 0 0 0;   
    load 186 0 [expr  3*$P3/32] 0 0 0 0;  
    load 187 0 [expr  3*$P3/32] 0 0 0 0;   
    load 188 0 [expr  5*$P3/48] 0 0 0 0; 
    load 189 0 [expr  5*$P3/48] 0 0 0 0;  
    load 190 0 [expr  5*$P3/48] 0 0 0 0;   
    load 191 0 [expr  3*$P3/32] 0 0 0 0;  
    load 192 0 [expr  3*$P3/32] 0 0 0 0; 
    load 193 [expr $P3/6] 0 0 0 0 0;   
    load 194 [expr $P3/6] 0 0 0 0 0; 
    load 195 [expr $P3/6] 0 0 0 0 0;   
    load 196 [expr $P3/6] 0 0 0 0 0; 
    load 197 [expr $P3/6] 0 0 0 0 0;   
    load 198 [expr $P3/6] 0 0 0 0 0;   
  }


# -------------------------------------------------------------------------
# Create ModelBuilder for 3D elements in X-Z Plane (with two-dimensions and 2 DOF/node)
# -------------------------------------------------------------------------

  model basic -ndm 2 -ndf 2  

# Create nodes & add to Domain - command: node nodeId xCrd yCrd zCrd
# Nodes 1 < node < 66 are for plane stress elements

  node  1 0 0;
  node  2 $deltaAB 0;
  node  7 $x1 0;
  node  8 $x2 0;
  node  9 $x3 0;
  node 11 0 $cH1;
  node 12 $deltaAB $cH1;
  node 17 $x1 $cH1;
  node 18 $x2 $cH1;
  node 19 $x3 $cH1;
  node 21 0 $cH2;
  node 22 $deltaAB $cH2;
  node 27 $x1 $cH2;
  node 28 $x2 $cH2;
  node 29 $x3 $cH2;
  node 31 0 $cH3;
  node 32 $deltaAB $cH3;
  node 37 $x1 $cH3;
  node 38 $x2 $cH3;
  node 39 $x3 $cH3;
  node 41 0 $cH4;
  node 42 $deltaAB $cH4;
  node 47 $x1 $cH4;
  node 48 $x2 $cH4;
  node 49 $x3 $cH4;
  node 51 0 $cH5;
  node 52 $deltaAB $cH5;
  node 57 $x1 $cH5;
  node 58 $x2 $cH5;
  node 59 $x3 $cH5;


# Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt? zRestrnt?

  fix 1 1 1;  
  fix 2 1 1; 


# tie nodes between beam, column, and wall elements

  equalDOF 11  111  1 3;
  equalDOF 12  112  1 3;
  equalDOF 21  121  1 3;
  equalDOF 22  122  1 3;
  equalDOF 27  167  1 3;
  equalDOF 28  168  1 3;
  equalDOF 29  169  1 3;
  equalDOF 31  131  1 3;
  equalDOF 32  132  1 3;
  equalDOF 41  141  1 3;
  equalDOF 42  142  1 3;
  equalDOF 51  151  1 3;
  equalDOF 52  152  1 3;
  equalDOF 57  183  1 3;
  equalDOF 58  184  1 3;
  equalDOF 59  185  1 3;

# -----------------------------------------------------------------
#  Define materials for ReinforceConcretePlaneStress element
# -----------------------------------------------------------------


# set fc fy E
  set wfc $fcw;
  set wfy $fy;
  set wE  200000.0;
  set rouw 0.0015;   # wall reinforcement ratio
  set ec 0.0025;

# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
  uniaxialMaterial    SteelZ01  6   $wfy    $wE  $wfc  $rouw


# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#  uniaxialMaterial ConcreteZ01  7 [expr -$wfc] -0.0025  
#  uniaxialMaterial ConcreteZ01  8 [expr -$wfc] -0.0025 

  uniaxialMaterial ConcreteZ02  7 [expr -$wfc] -0.0025  15000
  uniaxialMaterial ConcreteZ02  8 [expr -$wfc] -0.0025  15000


  set pi 3.141592654;
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2        rou1  rou2     fpc  fy  E0
  nDMaterial ReinforceConcretePlaneStress  9  0.0  6  6  7  7 [expr 0.0*$pi] [expr 0.5*$pi]  $rouw $rouw  $wfc $wfy $wE $ec


# -----------------------------------------------------------------
#  Define ReinforceConcretePlaneStress element
# -----------------------------------------------------------------

  set wtcc 150;    # Wall thickness C-C
  set wt11 100;      # wall thickness 1-1  

# Elements 2000 < element < 3000 are used for walls

  set k 0;
  while {$k < 5} {    
  # element quad eleID            node1          node2          node3           node4               thick     type         matID
    element quad [expr 10*$k+2001] [expr 10*$k+1] [expr 10*$k+7] [expr 10*$k+17] [expr 10*$k+11]     $wt11     PlaneStress  9
    element quad [expr 10*$k+2007] [expr 10*$k+7] [expr 10*$k+8] [expr 10*$k+18] [expr 10*$k+17]     $wt11     PlaneStress  9
    element quad [expr 10*$k+2008] [expr 10*$k+8] [expr 10*$k+9] [expr 10*$k+19] [expr 10*$k+18]     $wt11     PlaneStress  9 
    element quad [expr 10*$k+2009] [expr 10*$k+9] [expr 10*$k+2] [expr 10*$k+12] [expr 10*$k+19]     $wt11     PlaneStress  9 
   set k [expr $k+1]
  }

# -------------------------------------------------------------------------
# Create ModelBuilder for 3D elements in Y-Z Plane (with two-dimensions and 2 DOF/node)
# -------------------------------------------------------------------------

  model basic -ndm 2 -ndf 2  

# Create nodes & add to Domain - command: node nodeId xCrd yCrd zCrd
# Nodes 1 < node < 66 are for plane stress elements

  node  3  0 0;
  node  6  $delta12 0;
  node 13  0 $cH1;
  node 16  $delta12 $cH1;
  node 23  0 $cH2;
  node 26  $delta12 $cH2;
  node 33  0 $cH3;
  node 36  $delta12 $cH3;
  node 43  0 $cH4;
  node 46  $delta12 $cH4;
  node 53  0 $cH5;
  node 56  $delta12 $cH5;


# Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt? zRestrnt?

  fix 3 1 1;  
  fix 6 1 1; 


# tie nodes between beam, column, and wall elements

  equalDOF 13  113  2 3;
  equalDOF 16  116  2 3;
  equalDOF 23  123  2 3;
  equalDOF 26  126  2 3;
  equalDOF 33  133  2 3;
  equalDOF 36  136  2 3;
  equalDOF 43  143  2 3;
  equalDOF 46  146  2 3;
  equalDOF 53  153  2 3;
  equalDOF 56  156  2 3;
  
# -----------------------------------------------------------------
#  Define materials for ReinforceConcretePlaneStress element
# -----------------------------------------------------------------


# set fc fy E
  set wfc $fcw;
  set wfy $fy;
  set wE  200000.0;
  set rouw 0.0015;   # wall reinforcement ratio
  set ec 0.0025;

# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
  uniaxialMaterial    SteelZ01  6   $wfy    $wE  $wfc  $rouw


# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#  uniaxialMaterial ConcreteZ01  7 [expr -$wfc] -0.0025  
#  uniaxialMaterial ConcreteZ01  8 [expr -$wfc] -0.0025 

  uniaxialMaterial ConcreteZ02  7 [expr -$wfc] -0.0025  15000
  uniaxialMaterial ConcreteZ02  8 [expr -$wfc] -0.0025  15000


  set pi 3.141592654;
# NDMaterial: ReinforceConcretePlaneStress
#                                        tag  rho s1 s2 c1 c2    angle1         angle2        rou1  rou2     fpc  fy  E0
  nDMaterial ReinforceConcretePlaneStress  9  0.0  6  6  7  7 [expr 0.0*$pi] [expr 0.5*$pi]  $rouw $rouw  $wfc $wfy $wE $ec

# -----------------------------------------------------------------
#  Define ReinforceConcretePlaneStress element
# -----------------------------------------------------------------

# Elements 2000 < element < 3000 are used for walls


  set l 0;
  while {$l < 5} {    
  # element quad eleID            node1          node2          node3           node4               thick     type         matID
    element quad [expr 10*$l+2003] [expr 10*$l+3] [expr 10*$l+6] [expr 10*$l+16] [expr 10*$l+13]     $wt11     PlaneStress  9
   set l [expr $l+1]
  } 


# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------
# Start of analysis generation
# ------------------------------


# Create the system of equation, a sparse solver with partial pivoting
  system BandGeneral


# Create the constraint handler
  constraints Plain

# Create the DOF numberer
  numberer Plain


# Create the convergence test
  test NormDispIncrVaryIter 0.01 16 5 numStep 100 400 300 400 700 1020 1320 1540 1820 2100 2300 2570 2770 3050 1500 10 numIter 100 0 0 0 0 0 0 0 0 0 0 0 100 100 100 0  

# Create the solution algorithm
#   algorithm Newton 
#   algorithm NewtonLineSearch 0.8
    algorithm KrylovNewton

# Create the integration scheme, the DisplacementControl scheme
  integrator DisplacementPath 192 2 16 numStep 100 400 300 400 700 1020 1320 1540 1820 2100 2300 2570 2770 3050 1500 10 increment 0.001 -0.005 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01 0.01 -0.01


# Create the analysis object
  analysis Static


# initialize in case we need to do an initial stiffness iteration
  initialize

# ------------------------------
# End of analysis generation
# ------------------------------

# Create a recorder to monitor nodal displacements
  recorder Node -file NCREE_3D_Frame.out -time -node 192 -dof 2 disp

# perform the analysis
    analyze 21890
#    analyze 1500

# Print out the state of nodes
  print node 192

# Print out the state of elements
#  print ele 5 6 7 8