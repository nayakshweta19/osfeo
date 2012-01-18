# ----------------------------
# Start of model generation
# ----------------------------
wipe;
source Units&Constants.tcl;
# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 3 -ndf 6

# create the material
nDMaterial ElasticIsotropic   1   4.0e4   0.2  6.75e-6;

# -----------------------------------------------------------------
#  Define materials for 2D ReinforceConcretePlaneStress element
# -----------------------------------------------------------------
 
# set fc fy E
set wfc 34.5;  #[expr (6880.+6600.+7700.+6850.)/4.0*$psi/$MPa];
set wfy 414.0; #[expr (76.+63.+60)/3.0*$ksi/$MPa];
set wE  1.95e5; #[expr (29000+28000+28700)/3.0*$ksi/$MPa];
set rou1 0.0259;
set rou2 0.0226;
 
# UniaxialMaterial: steelZ01
#                               tag   fy     E0  fpc     rou
uniaxialMaterial    SteelZ01  11   $wfy    $wE  $wfc  $rou1
uniaxialMaterial    SteelZ01  12   $wfy    $wE  $wfc  $rou2
#uniaxialMaterial    SteelZ02  11   $wfy    $wE  $wfc  $rou1
#uniaxialMaterial    SteelZ02  12   $wfy    $wE  $wfc  $rou2


# UniaxialMaterial: concreteZ01
# ConcreteZ01                tag   f'c     ec0   
#uniaxialMaterial ConcreteZ01  13 [expr -$wfc] -0.0025 
#uniaxialMaterial ConcreteZ01  14 [expr -$wfc] -0.0025 
uniaxialMaterial ConcreteL02  13 [expr -$wfc] -0.0025 
uniaxialMaterial ConcreteL02  14 [expr -$wfc] -0.0025 
#uniaxialMaterial Concrete09 13 -20.1  -0.0021  -5.6  -0.03  0.14  +2.6  +1300 
#uniaxialMaterial Concrete09 14 -20.1  -0.0021  -5.6  -0.03  0.14  +2.6  +1300 

set pi 3.141592654;
# NDMaterial: ReinforceConcretePlateFiber
#                                        tag    rho   s1 s2 c1 c2    angle1         angle2         rou1   rou2  fpc  fy  E0
#nDMaterial FAReinforcedConcretePlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

#nDMaterial RAReinforcedConcretePlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002
nDMaterial CSMMRCPlateFiber 2  2.65e-6 11 12 13 14 [expr 0.0*$pi]  [expr 0.5*$pi]   $rou1 $rou2  $wfc $wfy $wE 0.002

# ---------------------------------------
# Define Section tag
# ---------------------------------------
#section PlateFiber secTag  ndMatTag  h
set elastEleSecTag 2;
set thick [expr $in*6.0/$mm];
section  PlateFiber  $elastEleSecTag   1  $thick

set eleSecTag 3;
section  PlateFiber  $eleSecTag   2  $thick  

# ---------------------------------------
# Define geometry
# ---------------------------------------

# define some  parameters

set mass 0.0;

set xWidth [expr 7.8E+1*$in/$mm];
set yWidth [expr 5.4E+1*$in/$mm];
set height [expr 2.88E+2*$in/$mm];

set numX 4; # num elements in x direction
set numY 4; # num elements in y direction ,divied by 2, for middle wall
set numZ 16; # num elements in z direction ,divied by 4, for 4 story\
#     ^
#     | Y
#       
#     |
#     |
#     |
#     |--------------  -> X  
#     |
#     |
#     |
#
set nodeNum 1

set incrX [expr $xWidth/(1.0*$numX)]
set incrY [expr $yWidth/(1.0*$numY)]
set incrZ [expr $height/(1.0*$numZ)]
#web wall

node    1          1981.2               0          7315.2
node    2          1981.2               0          7086.6
node    3         1733.55               0          7315.2
node    4         1733.55               0          7086.6
node    5          1981.2               0            6858
node    6          1485.9               0          7315.2
node    7         1733.55               0            6858
node    8          1485.9               0          7086.6
node    9          1485.9               0            6858
node   10          1981.2               0          6629.4
node   11         1733.55               0          6629.4
node   12         1238.25               0          7315.2
node   13         1238.25               0          7086.6
node   14          1485.9               0          6629.4
node   15         1238.25               0            6858
node   16          1981.2               0          6400.8
node   17         1733.55               0          6400.8
node   18           990.6               0          7315.2
node   19         1238.25               0          6629.4
node   20           990.6               0          7086.6
node   21          1485.9               0          6400.8
node   22           990.6               0            6858
node   23          1981.2               0          6172.2
node   24         1733.55               0          6172.2
node   25         1238.25               0          6400.8
node   26           990.6               0          6629.4
node   27          742.95               0          7315.2
node   28          1485.9               0          6172.2
node   29          742.95               0          7086.6
node   30          742.95               0            6858
node   31           990.6               0          6400.8
node   32         1238.25               0          6172.2
node   33          1981.2               0          5943.6
node   34         1733.55               0          5943.6
node   35          742.95               0          6629.4
node   36          1485.9               0          5943.6
node   37           495.3               0          7315.2
node   38           495.3               0          7086.6
node   39           990.6               0          6172.2
node   40          742.95               0          6400.8
node   41           495.3               0            6858
node   42         1238.25               0          5943.6
node   43          1981.2               0            5715
node   44         1733.55               0            5715
node   45           495.3               0          6629.4
node   46          1485.9               0            5715
node   47          742.95               0          6172.2
node   48           990.6               0          5943.6
node   49          247.65               0          7315.2
node   50           495.3               0          6400.8
node   51          247.65               0          7086.6
node   52         1238.25               0            5715
node   53          247.65               0            6858
node   54          1981.2               0          5486.4
node   55         1733.55               0          5486.4
node   56          742.95               0          5943.6
node   57          247.65               0          6629.4
node   58           495.3               0          6172.2
node   59           990.6               0            5715
node   60          1485.9               0          5486.4
node   61          247.65               0          6400.8
node   62         1238.25               0          5486.4
node   63               0               0          7315.2
node   64               0          171.45          7315.2
node   65               0         -171.45          7315.2
node   66               0               0          7086.6
node   67               0         -171.45          7086.6
node   68               0          171.45          7086.6
node   69               0           342.9          7315.2
node   70               0          -342.9          7315.2
node   71           495.3               0          5943.6
node   72          742.95               0            5715
node   73               0           342.9          7086.6
node   74               0          -342.9          7086.6
node   75               0               0            6858
node   76               0         -171.45            6858
node   77               0          171.45            6858
node   78               0          514.35          7315.2
node   79               0         -514.35          7315.2
node   80          1981.2               0          5257.8
node   81               0          514.35          7086.6
node   82               0         -514.35          7086.6
node   83               0           342.9            6858
node   84               0          -342.9            6858
node   85         1733.55               0          5257.8
node   86          247.65               0          6172.2
node   87           990.6               0          5486.4
node   88               0          -685.8          7315.2
node   89               0           685.8          7315.2
node   90               0               0          6629.4
node   91               0          514.35            6858
node   92               0         -514.35            6858
node   93               0          171.45          6629.4
node   94               0         -171.45          6629.4
node   95               0           685.8          7086.6
node   96               0          -685.8          7086.6
node   97          1485.9               0          5257.8
node   98               0           342.9          6629.4
node   99               0          -342.9          6629.4
node  100               0           685.8            6858
node  101               0          -685.8            6858
node  102               0         -514.35          6629.4
node  103               0          514.35          6629.4
node  104               0               0          6400.8
node  105           495.3               0            5715
node  106         1238.25               0          5257.8
node  107               0         -171.45          6400.8
node  108               0          171.45          6400.8
node  109               0           685.8          6629.4
node  110               0          -685.8          6629.4
node  111          742.95               0          5486.4
node  112               0          -342.9          6400.8
node  113               0           342.9          6400.8
node  114          247.65               0          5943.6
node  115               0         -514.35          6400.8
node  116               0          514.35          6400.8
node  117           990.6               0          5257.8
node  118          1981.2               0          5029.2
node  119               0               0          6172.2
node  120               0           685.8          6400.8
node  121               0          -685.8          6400.8
node  122               0         -171.45          6172.2
node  123               0          171.45          6172.2
node  124         1733.55               0          5029.2
node  125               0           342.9          6172.2
node  126               0          -342.9          6172.2
node  127          1485.9               0          5029.2
node  128               0         -514.35          6172.2
node  129               0          514.35          6172.2
node  130           495.3               0          5486.4
node  131          247.65               0            5715
node  132               0          -685.8          6172.2
node  133               0           685.8          6172.2
node  134          742.95               0          5257.8
node  135         1238.25               0          5029.2
node  136               0               0          5943.6
node  137               0         -171.45          5943.6
node  138               0          171.45          5943.6
node  139               0          -342.9          5943.6
node  140               0           342.9          5943.6
node  141               0          514.35          5943.6
node  142               0         -514.35          5943.6
node  143           990.6               0          5029.2
node  144               0          -685.8          5943.6
node  145               0           685.8          5943.6
node  146          1981.2               0          4800.6
node  147          247.65               0          5486.4
node  148         1733.55               0          4800.6
node  149           495.3               0          5257.8
node  150               0               0            5715
node  151               0         -171.45            5715
node  152               0          171.45            5715
node  153          1485.9               0          4800.6
node  154               0          -342.9            5715
node  155               0           342.9            5715
node  156               0          514.35            5715
node  157               0         -514.35            5715
node  158          742.95               0          5029.2
node  159         1238.25               0          4800.6
node  160               0          -685.8            5715
node  161               0           685.8            5715
node  162          247.65               0          5257.8
node  163               0               0          5486.4
node  164               0          171.45          5486.4
node  165               0         -171.45          5486.4
node  166           990.6               0          4800.6
node  167               0          -342.9          5486.4
node  168               0           342.9          5486.4
node  169           495.3               0          5029.2
node  170          1981.2               0            4572
node  171               0         -514.35          5486.4
node  172               0          514.35          5486.4
node  173         1733.55               0            4572
node  174               0           685.8          5486.4
node  175               0          -685.8          5486.4
node  176          1485.9               0            4572
node  177          742.95               0          4800.6
node  178         1238.25               0            4572
node  179               0               0          5257.8
node  180               0         -171.45          5257.8
node  181               0          171.45          5257.8
node  182          247.65               0          5029.2
node  183               0          -342.9          5257.8
node  184               0           342.9          5257.8
node  185               0         -514.35          5257.8
node  186               0          514.35          5257.8
node  187           990.6               0            4572
node  188           495.3               0          4800.6
node  189               0          -685.8          5257.8
node  190               0           685.8          5257.8
node  191          1981.2               0          4343.4
node  192         1733.55               0          4343.4
node  193          742.95               0            4572
node  194          1485.9               0          4343.4
node  195               0               0          5029.2
node  196               0          171.45          5029.2
node  197               0         -171.45          5029.2
node  198               0           342.9          5029.2
node  199               0          -342.9          5029.2
node  200          247.65               0          4800.6
node  201         1238.25               0          4343.4
node  202               0          514.35          5029.2
node  203               0         -514.35          5029.2
node  204               0           685.8          5029.2
node  205               0          -685.8          5029.2
node  206           495.3               0            4572
node  207           990.6               0          4343.4
node  208          1981.2               0          4114.8
node  209               0               0          4800.6
node  210               0          171.45          4800.6
node  211               0         -171.45          4800.6
node  212         1733.55               0          4114.8
node  213          742.95               0          4343.4
node  214               0          -342.9          4800.6
node  215               0           342.9          4800.6
node  216          1485.9               0          4114.8
node  217               0         -514.35          4800.6
node  218               0          514.35          4800.6
node  219          247.65               0            4572
node  220               0           685.8          4800.6
node  221               0          -685.8          4800.6
node  222         1238.25               0          4114.8
node  223           495.3               0          4343.4
node  224           990.6               0          4114.8
node  225               0               0            4572
node  226               0         -171.45            4572
node  227               0          171.45            4572
node  228               0           342.9            4572
node  229               0          -342.9            4572
node  230               0         -514.35            4572
node  231               0          514.35            4572
node  232          1981.2               0          3886.2
node  233          742.95               0          4114.8
node  234         1733.55               0          3886.2
node  235          247.65               0          4343.4
node  236               0          -685.8            4572
node  237               0           685.8            4572
node  238          1485.9               0          3886.2
node  239         1238.25               0          3886.2
node  240           495.3               0          4114.8
node  241           990.6               0          3886.2
node  242               0               0          4343.4
node  243               0         -171.45          4343.4
node  244               0          171.45          4343.4
node  245               0          -342.9          4343.4
node  246               0           342.9          4343.4
node  247               0         -514.35          4343.4
node  248               0          514.35          4343.4
node  249               0          -685.8          4343.4
node  250               0           685.8          4343.4
node  251          247.65               0          4114.8
node  252          742.95               0          3886.2
node  253          1981.2               0          3657.6
node  254         1733.55               0          3657.6
node  255          1485.9               0          3657.6
node  256         1238.25               0          3657.6
node  257           495.3               0          3886.2
node  258               0               0          4114.8
node  259               0          171.45          4114.8
node  260               0         -171.45          4114.8
node  261               0           342.9          4114.8
node  262               0          -342.9          4114.8
node  263           990.6               0          3657.6
node  264               0          514.35          4114.8
node  265               0         -514.35          4114.8
node  266               0           685.8          4114.8
node  267               0          -685.8          4114.8
node  268          247.65               0          3886.2
node  269          742.95               0          3657.6
node  270          1981.2               0            3429
node  271         1733.55               0            3429
node  272          1485.9               0            3429
node  273           495.3               0          3657.6
node  274         1238.25               0            3429
node  275               0               0          3886.2
node  276               0          171.45          3886.2
node  277               0         -171.45          3886.2
node  278               0          -342.9          3886.2
node  279               0           342.9          3886.2
node  280               0         -514.35          3886.2
node  281               0          514.35          3886.2
node  282           990.6               0            3429
node  283               0          -685.8          3886.2
node  284               0           685.8          3886.2
node  285          247.65               0          3657.6
node  286          742.95               0            3429
node  287          1981.2               0          3200.4
node  288         1733.55               0          3200.4
node  289          1485.9               0          3200.4
node  290               0               0          3657.6
node  291           495.3               0            3429
node  292               0          171.45          3657.6
node  293               0         -171.45          3657.6
node  294               0          -342.9          3657.6
node  295               0           342.9          3657.6
node  296         1238.25               0          3200.4
node  297               0         -514.35          3657.6
node  298               0          514.35          3657.6
node  299               0          -685.8          3657.6
node  300               0           685.8          3657.6
node  301           990.6               0          3200.4
node  302          247.65               0            3429
node  303          742.95               0          3200.4
node  304          1981.2               0          2971.8
node  305         1733.55               0          2971.8
node  306               0               0            3429
node  307               0          171.45            3429
node  308               0         -171.45            3429
node  309          1485.9               0          2971.8
node  310           495.3               0          3200.4
node  311               0           342.9            3429
node  312               0          -342.9            3429
node  313               0         -514.35            3429
node  314               0          514.35            3429
node  315         1238.25               0          2971.8
node  316               0          -685.8            3429
node  317               0           685.8            3429
node  318           990.6               0          2971.8
node  319          247.65               0          3200.4
node  320          742.95               0          2971.8
node  321               0               0          3200.4
node  322               0         -171.45          3200.4
node  323               0          171.45          3200.4
node  324          1981.2               0          2743.2
node  325         1733.55               0          2743.2
node  326               0          -342.9          3200.4
node  327               0           342.9          3200.4
node  328           495.3               0          2971.8
node  329               0         -514.35          3200.4
node  330               0          514.35          3200.4
node  331          1485.9               0          2743.2
node  332               0          -685.8          3200.4
node  333               0           685.8          3200.4
node  334         1238.25               0          2743.2
node  335          247.65               0          2971.8
node  336           990.6               0          2743.2
node  337          742.95               0          2743.2
node  338               0               0          2971.8
node  339               0         -171.45          2971.8
node  340               0          171.45          2971.8
node  341               0           342.9          2971.8
node  342               0          -342.9          2971.8
node  343          1981.2               0          2514.6
node  344               0         -514.35          2971.8
node  345               0          514.35          2971.8
node  346         1733.55               0          2514.6
node  347           495.3               0          2743.2
node  348               0          -685.8          2971.8
node  349               0           685.8          2971.8
node  350          1485.9               0          2514.6
node  351         1238.25               0          2514.6
node  352          247.65               0          2743.2
node  353           990.6               0          2514.6
node  354          742.95               0          2514.6
node  355               0               0          2743.2
node  356               0          171.45          2743.2
node  357               0         -171.45          2743.2
node  358               0           342.9          2743.2
node  359               0          -342.9          2743.2
node  360               0         -514.35          2743.2
node  361               0          514.35          2743.2
node  362           495.3               0          2514.6
node  363          1981.2               0            2286
node  364               0          -685.8          2743.2
node  365               0           685.8          2743.2
node  366         1733.55               0            2286
node  367          1485.9               0            2286
node  368         1238.25               0            2286
node  369          247.65               0          2514.6
node  370           990.6               0            2286
node  371          742.95               0            2286
node  372               0               0          2514.6
node  373               0          171.45          2514.6
node  374               0         -171.45          2514.6
node  375               0          -342.9          2514.6
node  376               0           342.9          2514.6
node  377               0         -514.35          2514.6
node  378               0          514.35          2514.6
node  379               0           685.8          2514.6
node  380               0          -685.8          2514.6
node  381           495.3               0            2286
node  382          1981.2               0          2057.4
node  383         1733.55               0          2057.4
node  384          1485.9               0          2057.4
node  385         1238.25               0          2057.4
node  386          247.65               0            2286
node  387           990.6               0          2057.4
node  388          742.95               0          2057.4
node  389               0               0            2286
node  390               0         -171.45            2286
node  391               0          171.45            2286
node  392               0           342.9            2286
node  393               0          -342.9            2286
node  394               0         -514.35            2286
node  395               0          514.35            2286
node  396               0           685.8            2286
node  397               0          -685.8            2286
node  398           495.3               0          2057.4
node  399          1981.2               0          1828.8
node  400         1733.55               0          1828.8
node  401          1485.9               0          1828.8
node  402          247.65               0          2057.4
node  403         1238.25               0          1828.8
node  404           990.6               0          1828.8
node  405               0               0          2057.4
node  406               0         -171.45          2057.4
node  407               0          171.45          2057.4
node  408          742.95               0          1828.8
node  409               0          -342.9          2057.4
node  410               0           342.9          2057.4
node  411               0         -514.35          2057.4
node  412               0          514.35          2057.4
node  413               0           685.8          2057.4
node  414               0          -685.8          2057.4
node  415           495.3               0          1828.8
node  416          1981.2               0          1600.2
node  417         1733.55               0          1600.2
node  418          1485.9               0          1600.2
node  419          247.65               0          1828.8
node  420         1238.25               0          1600.2
node  421           990.6               0          1600.2
node  422               0               0          1828.8
node  423               0         -171.45          1828.8
node  424               0          171.45          1828.8
node  425               0          -342.9          1828.8
node  426               0           342.9          1828.8
node  427          742.95               0          1600.2
node  428               0          514.35          1828.8
node  429               0         -514.35          1828.8
node  430               0          -685.8          1828.8
node  431               0           685.8          1828.8
node  432           495.3               0          1600.2
node  433          1981.2               0          1371.6
node  434         1733.55               0          1371.6
node  435          1485.9               0          1371.6
node  436          247.65               0          1600.2
node  437         1238.25               0          1371.6
node  438           990.6               0          1371.6
node  439               0               0          1600.2
node  440               0         -171.45          1600.2
node  441               0          171.45          1600.2
node  442               0           342.9          1600.2
node  443               0          -342.9          1600.2
node  444               0         -514.35          1600.2
node  445               0          514.35          1600.2
node  446          742.95               0          1371.6
node  447               0          -685.8          1600.2
node  448               0           685.8          1600.2
node  449           495.3               0          1371.6
node  450          1981.2               0            1143
node  451         1733.55               0            1143
node  452          247.65               0          1371.6
node  453          1485.9               0            1143
node  454         1238.25               0            1143
node  455           990.6               0            1143
node  456               0               0          1371.6
node  457               0          171.45          1371.6
node  458               0         -171.45          1371.6
node  459               0           342.9          1371.6
node  460               0          -342.9          1371.6
node  461               0         -514.35          1371.6
node  462               0          514.35          1371.6
node  463          742.95               0            1143
node  464               0          -685.8          1371.6
node  465               0           685.8          1371.6
node  466           495.3               0            1143
node  467          1981.2               0           914.4
node  468         1733.55               0           914.4
node  469          247.65               0            1143
node  470          1485.9               0           914.4
node  471         1238.25               0           914.4
node  472           990.6               0           914.4
node  473               0               0            1143
node  474               0          171.45            1143
node  475               0         -171.45            1143
node  476               0           342.9            1143
node  477               0          -342.9            1143
node  478               0         -514.35            1143
node  479               0          514.35            1143
node  480               0          -685.8            1143
node  481               0           685.8            1143
node  482          742.95               0           914.4
node  483           495.3               0           914.4
node  484          1981.2               0           685.8
node  485          247.65               0           914.4
node  486         1733.55               0           685.8
node  487          1485.9               0           685.8
node  488         1238.25               0           685.8
node  489               0               0           914.4
node  490               0         -171.45           914.4
node  491               0          171.45           914.4
node  492           990.6               0           685.8
node  493               0           342.9           914.4
node  494               0          -342.9           914.4
node  495               0          514.35           914.4
node  496               0         -514.35           914.4
node  497               0          -685.8           914.4
node  498               0           685.8           914.4
node  499          742.95               0           685.8
node  500           495.3               0           685.8
node  501          247.65               0           685.8
node  502          1981.2               0           457.2
node  503         1733.55               0           457.2
node  504          1485.9               0           457.2
node  505         1238.25               0           457.2
node  506               0               0           685.8
node  507               0         -171.45           685.8
node  508               0          171.45           685.8
node  509               0          -342.9           685.8
node  510               0           342.9           685.8
node  511           990.6               0           457.2
node  512               0         -514.35           685.8
node  513               0          514.35           685.8
node  514               0          -685.8           685.8
node  515               0           685.8           685.8
node  516          742.95               0           457.2
node  517           495.3               0           457.2
node  518          247.65               0           457.2
node  519          1981.2               0           228.6
node  520         1733.55               0           228.6
node  521          1485.9               0           228.6
node  522         1238.25               0           228.6
node  523               0               0           457.2
node  524               0          171.45           457.2
node  525               0         -171.45           457.2
node  526               0           342.9           457.2
node  527               0          -342.9           457.2
node  528           990.6               0           228.6
node  529               0          514.35           457.2
node  530               0         -514.35           457.2
node  531               0           685.8           457.2
node  532               0          -685.8           457.2
node  533          742.95               0           228.6
node  534           495.3               0           228.6
node  535          247.65               0           228.6
node  536          1981.2               0               0
node  537         1733.55               0               0
node  538          1485.9               0               0
node  539         1238.25               0               0
node  540               0               0           228.6
node  541               0          171.45           228.6
node  542               0         -171.45           228.6
node  543               0           342.9           228.6
node  544               0          -342.9           228.6
node  545               0          514.35           228.6
node  546               0         -514.35           228.6
node  547           990.6               0               0
node  548               0           685.8           228.6
node  549               0          -685.8           228.6
node  550          742.95               0               0
node  551           495.3               0               0
node  552          247.65               0               0
node  553               0               0               0
node  554               0          171.45               0
node  555               0         -171.45               0
node  556               0           342.9               0
node  557               0          -342.9               0
node  558               0          514.35               0
node  559               0         -514.35               0
node  560               0           685.8               0
node  561               0          -685.8               0

element shellNL 1 37 41 75 63 38 53 66 49 51		$eleSecTag
element shellNL 2 18 22 41 37 20 30 38 27 29		$eleSecTag
element shellNL 3 6 9 22 18 8 15 20 12 13		$eleSecTag
element shellNL 4 1 5 9 6 2 7 8 3 4			$eleSecTag
element shellNL 5 41 50 104 75 45 61 90 53 57		$eleSecTag
element shellNL 6 22 31 50 41 26 40 45 30 35		$eleSecTag
element shellNL 7 9 21 31 22 14 25 26 15 19		$eleSecTag
element shellNL 8 5 16 21 9 10 17 14 7 11		$eleSecTag
element shellNL 9 50 71 136 104 58 114 119 61 86	$eleSecTag
element shellNL 10 31 48 71 50 39 56 58 40 47		$eleSecTag
element shellNL 11 21 36 48 31 28 42 39 25 32		$eleSecTag
element shellNL 12 16 33 36 21 23 34 28 17 24		$eleSecTag
element shellNL 13 71 130 163 136 105 147 150 114 131	$eleSecTag
element shellNL 14 48 87 130 71 59 111 105 56 72	$eleSecTag
element shellNL 15 36 60 87 48 46 62 59 42 52		$eleSecTag
element shellNL 16 33 54 60 36 43 55 46 34 44		$eleSecTag
element shellNL 17 551 517 523 553 534 518 540 552 535	$eleSecTag
element shellNL 18 547 511 517 551 528 516 534 550 533	$eleSecTag
element shellNL 19 538 504 511 547 521 505 528 539 522	$eleSecTag
element shellNL 20 536 502 504 538 519 503 521 537 520	$eleSecTag
element shellNL 21 517 483 489 523 500 485 506 518 501	$eleSecTag
element shellNL 22 511 472 483 517 492 482 500 516 499	$eleSecTag
element shellNL 23 504 470 472 511 487 471 492 505 488	$eleSecTag
element shellNL 24 502 467 470 504 484 468 487 503 486	$eleSecTag
element shellNL 25 483 449 456 489 466 452 473 485 469	$eleSecTag
element shellNL 26 472 438 449 483 455 446 466 482 463	$eleSecTag
element shellNL 27 470 435 438 472 453 437 455 471 454	$eleSecTag
element shellNL 28 467 433 435 470 450 434 453 468 451	$eleSecTag
element shellNL 29 449 415 422 456 432 419 439 452 436	$eleSecTag
element shellNL 30 438 404 415 449 421 408 432 446 427	$eleSecTag
element shellNL 31 435 401 404 438 418 403 421 437 420	$eleSecTag
element shellNL 32 433 399 401 435 416 400 418 434 417	$eleSecTag
element shellNL 33 401 367 363 399 384 366 382 400 383	$eleSecTag
element shellNL 34 404 370 367 401 387 368 384 403 385	$eleSecTag
element shellNL 35 415 381 370 404 398 371 387 408 388	$eleSecTag
element shellNL 36 422 389 381 415 405 386 398 419 402	$eleSecTag
element shellNL 37 367 331 324 363 350 325 343 366 346	$eleSecTag
element shellNL 38 370 336 331 367 353 334 350 368 351	$eleSecTag
element shellNL 39 381 347 336 370 362 337 353 371 354	$eleSecTag
element shellNL 40 389 355 347 381 372 352 362 386 369	$eleSecTag
element shellNL 41 331 289 287 324 309 288 304 325 305	$eleSecTag
element shellNL 42 336 301 289 331 318 296 309 334 315	$eleSecTag
element shellNL 43 347 310 301 336 328 303 318 337 320	$eleSecTag
element shellNL 44 355 321 310 347 338 319 328 352 335	$eleSecTag
element shellNL 45 289 255 253 287 272 254 270 288 271	$eleSecTag
element shellNL 46 301 263 255 289 282 256 272 296 274	$eleSecTag
element shellNL 47 310 273 263 301 291 269 282 303 286	$eleSecTag
element shellNL 48 321 290 273 310 306 285 291 319 302	$eleSecTag
element shellNL 49 456 459 426 422 457 442 424 439 441	$eleSecTag
element shellNL 50 489 493 459 456 491 476 457 473 474	$eleSecTag
element shellNL 51 523 526 493 489 524 510 491 506 508	$eleSecTag
element shellNL 52 553 556 526 523 554 543 524 540 541	$eleSecTag
element shellNL 53 459 465 431 426 462 448 428 442 445	$eleSecTag
element shellNL 54 493 498 465 459 495 481 462 476 479	$eleSecTag
element shellNL 55 526 531 498 493 529 515 495 510 513	$eleSecTag
element shellNL 56 556 560 531 526 558 548 529 543 545	$eleSecTag
element shellNL 57 456 460 425 422 458 443 423 439 440	$eleSecTag
element shellNL 58 489 494 460 456 490 477 458 473 475	$eleSecTag
element shellNL 59 523 527 494 489 525 509 490 506 507	$eleSecTag
element shellNL 60 553 557 527 523 555 544 525 540 542	$eleSecTag
element shellNL 61 460 464 430 425 461 447 429 443 444	$eleSecTag
element shellNL 62 494 497 464 460 496 480 461 477 478	$eleSecTag
element shellNL 63 527 532 497 494 530 514 496 509 512	$eleSecTag
element shellNL 64 557 561 532 527 559 549 530 544 546	$eleSecTag
element shellNL 65 389 392 426 422 391 410 424 405 407	$eleSecTag
element shellNL 66 355 358 392 389 356 376 391 372 373	$eleSecTag
element shellNL 67 321 327 358 355 323 341 356 338 340	$eleSecTag
element shellNL 68 290 295 327 321 292 311 323 306 307	$eleSecTag
element shellNL 69 392 396 431 426 395 413 428 410 412	$eleSecTag
element shellNL 70 358 365 396 392 361 379 395 376 378	$eleSecTag
element shellNL 71 327 333 365 358 330 349 361 341 345	$eleSecTag
element shellNL 72 295 300 333 327 298 317 330 311 314	$eleSecTag
element shellNL 73 397 393 425 430 394 409 429 414 411	$eleSecTag
element shellNL 74 364 359 393 397 360 375 394 380 377	$eleSecTag
element shellNL 75 332 326 359 364 329 342 360 348 344	$eleSecTag
element shellNL 76 299 294 326 332 297 312 329 316 313	$eleSecTag
element shellNL 77 393 389 422 425 390 405 423 409 406	$eleSecTag
element shellNL 78 359 355 389 393 357 372 390 375 374	$eleSecTag
element shellNL 79 326 321 355 359 322 338 357 342 339	$eleSecTag
element shellNL 80 294 290 321 326 293 306 322 312 308	$eleSecTag
element shellNL 81 130 169 195 163 149 182 179 147 162	$eleSecTag
element shellNL 82 87 143 169 130 117 158 149 111 134	$eleSecTag
element shellNL 83 60 127 143 87 97 135 117 62 106	$eleSecTag
element shellNL 84 54 118 127 60 80 124 97 55 85	$eleSecTag
element shellNL 85 169 206 225 195 188 219 209 182 200	$eleSecTag
element shellNL 86 143 187 206 169 166 193 188 158 177	$eleSecTag
element shellNL 87 127 176 187 143 153 178 166 135 159	$eleSecTag
element shellNL 88 118 170 176 127 146 173 153 124 148	$eleSecTag
element shellNL 89 206 240 258 225 223 251 242 219 235	$eleSecTag
element shellNL 90 187 224 240 206 207 233 223 193 213	$eleSecTag
element shellNL 91 176 216 224 187 194 222 207 178 201	$eleSecTag
element shellNL 92 170 208 216 176 191 212 194 173 192	$eleSecTag
element shellNL 93 240 273 290 258 257 285 275 251 268	$eleSecTag
element shellNL 94 224 263 273 240 241 269 257 233 252	$eleSecTag
element shellNL 95 216 255 263 224 238 256 241 222 239	$eleSecTag
element shellNL 96 208 253 255 216 232 254 238 212 234	$eleSecTag
element shellNL 97 258 262 294 290 260 278 293 275 277	$eleSecTag
element shellNL 98 225 229 262 258 226 245 260 242 243	$eleSecTag
element shellNL 99 195 199 229 225 197 214 226 209 211	$eleSecTag
element shellNL 100 163 167 199 195 165 183 197 179 180	$eleSecTag
element shellNL 101 262 267 299 294 265 283 297 278 280	$eleSecTag
element shellNL 102 229 236 267 262 230 249 265 245 247	$eleSecTag
element shellNL 103 199 205 236 229 203 221 230 214 217	$eleSecTag
element shellNL 104 167 175 205 199 171 189 203 183 185	$eleSecTag
element shellNL 105 204 198 168 174 202 184 172 190 186	$eleSecTag
element shellNL 106 237 228 198 204 231 215 202 220 218	$eleSecTag
element shellNL 107 266 261 228 237 264 246 231 250 248	$eleSecTag
element shellNL 108 300 295 261 266 298 279 264 284 281	$eleSecTag
element shellNL 109 198 195 163 168 196 179 164 184 181	$eleSecTag
element shellNL 110 228 225 195 198 227 209 196 215 210	$eleSecTag
element shellNL 111 261 258 225 228 259 242 227 246 244	$eleSecTag
element shellNL 112 295 290 258 261 292 275 259 279 276	$eleSecTag
element shellNL 113 136 140 168 163 138 155 164 150 152	$eleSecTag
element shellNL 114 104 113 140 136 108 125 138 119 123	$eleSecTag
element shellNL 115 75 83 113 104 77 98 108 90 93	$eleSecTag
element shellNL 116 63 69 83 75 64 73 77 66 68		$eleSecTag
element shellNL 117 140 145 174 168 141 161 172 155 156	$eleSecTag
element shellNL 118 113 120 145 140 116 133 141 125 129	$eleSecTag
element shellNL 119 83 100 120 113 91 109 116 98 103	$eleSecTag
element shellNL 120 69 89 100 83 78 95 91 73 81		$eleSecTag
element shellNL 121 136 139 167 163 137 154 165 150 151	$eleSecTag
element shellNL 122 104 112 139 136 107 126 137 119 122	$eleSecTag
element shellNL 123 75 84 112 104 76 99 107 90 94	$eleSecTag
element shellNL 124 63 70 84 75 65 74 76 66 67		$eleSecTag
element shellNL 125 139 144 175 167 142 160 171 154 157	$eleSecTag
element shellNL 126 112 121 144 139 115 132 142 126 128	$eleSecTag
element shellNL 127 84 101 121 112 92 110 115 99 102	$eleSecTag
element shellNL 128 70 88 101 84 79 96 92 74 82		$eleSecTag

fixZ 0.0 1 1 1 1 1 1
# end web wall

#---------------------------------------------------
# Perform the eigenVector output
#---------------------------------------------------
set lambda [eigen 15];
record; puts "no graivity analysis period and frequency"
foreach lam $lambda {
  if { $lam > 0.0 } {
    lappend omega [expr sqrt($lam)]
    lappend f [expr sqrt($lam)/(2*$pi)]
    lappend T [expr (2*$pi)/sqrt($lam)]
    puts "Tn = [expr (2*$pi)/sqrt($lam)] s";
  } else {
    lappend omega -1.0
    lappend f "complex"
    lappend T "complex"
    puts "Tn = complex";
  }
}

timeSeries Linear 1
pattern  Plain       1    "Linear"  { 
puts $fileModel "pattern  Plain       1 {"
    # Load    nodeTag    LoadValues 
    for { set j 0 } { $j <= $numY} {incr j 1} {
      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 
      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 "
    }
    for { set j 1 } { $j <= $numX} {incr j 1} {
      load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 
      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000  +0.000000E+000  -3.000000E+004  +0.000000E+000  +0.000000E+000  +0.000000E+000 "
    }
    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
puts $fileModel "}"
}

# Constraint Handler Lagrange
constraints  Plain;  
# Convergence Test 
test  NormDispIncr  +1.000000E-003    250     0   
#test EnergyIncr  1.0e-12    10         0
# Integrator 
integrator  LoadControl  +1.000000E+000 
# Solution Algorithm 
algorithm  KrylovNewton 
# DOF Numberer 
numberer  Plain 
# System of Equations BandGeneral SparseGeneral UmfPackProfileSPD
system  BandGeneral 
# Analysis Type 
analysis  Static 
#---------------------------------------------------
# Perform the analysis
#---------------------------------------------------
set ok [analyze   10]
puts "Gravity Static Analysis OK= $ok (0 for ok, other for error)";

loadConst -time 0.0
#---------------------------------------------------
# Perform the eigenVector output
#---------------------------------------------------
set lambda [eigen 15];
record;
foreach lam $lambda {
  if { $lam > 0.0 } {
    lappend omega [expr sqrt($lam)]
    lappend f [expr sqrt($lam)/(2*$pi)]
    lappend T [expr (2*$pi)/sqrt($lam)]
    puts "Tn = [expr (2*$pi)/sqrt($lam)] s";
  } else {
    lappend omega -1.0
    lappend f "complex"
    lappend T "complex"
    puts "Tn = complex";
  }
}
set nodeIDs [getNodeTags];
set endNodeTag [lindex $nodeIDs end];
set eleIDs [getEleTags];
set endEleTag [lindex $eleIDs end];
# get values of eigenvectors for translational DOFs
set fileEig [open "eigFile.eig" "w"]
for { set k 1 } { $k <= 15 } { incr k } {
  set tempStr "[lindex $T [expr $k-1]]"
  foreach tg $nodeIDs {
    lappend tempStr [lindex [nodeEigenvector $tg $k 1] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 2] 0]
    lappend tempStr [lindex [nodeEigenvector $tg $k 3] 0]
  }
  puts $fileEig $tempStr;
}
close $fileEig 

# create the display
recorder display nodeNum 10 10 800 600 -wipe
# next three commmands define viewing system, all values in global coords
prp -80000 -10000 40000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 0 1;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -4000 4000 -4000 4000; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 -1 10; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# create the display
recorder display Deform 810 10 800 600 -wipe
# next three commmands define viewing system, all values in global coords
prp -10000 -100000 70000;          # eye location in local coord sys defined by viewing system
#vrp 0 -500 250;                  # point on the view plane in global coord, center of local viewing system
vup 0 0 1;                        # dirn defining up direction of view plane
#vpn -1 -1 0.5;                   # direction of outward normal to view plane
viewWindow -4000 4000 -4000 4000; # view bounds uMin, uMax, vMin, vMax in local coords
#plane 0 150;                     # distance to front and back clipping planes from eye
port -1 1 -1 1;                   # area of window that will be drawn into
projection 1;                     # projection mode
fill 0;                            # fill mode
display 1 5 1; # display -$nEigen 0 $dAmp;  # display mode shape for mode $nEigen
                 # display 1 -1 0 ;           # display node numbers
                 # display 1 5 $dAmp;         # display deformed shape

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
loadConst -time 0.0

# ---------------------------------- 
# Analysis: PushoverCase 
# ---------------------------------- 
# ---------------------------------- 
# Define time series 
# ---------------------------------- 

# TimeSeries "TimeSeriesX":    tsTag    dt    filePath    cFactor 
timeSeries  Path       2  -dt  +1.000000E-002  -filePath  TimeSeriesX.thf  -factor  +1.000000E-003 

# TimeSeries "TimeSeriesY":    tsTag    dt    filePath    cFactor 
timeSeries  Path       3  -dt  +1.000000E-002  -filePath  TimeSeriesY.thf  -factor  +1.000000E-003 

# LoadPattern "PushoverPattern":    patternTag    tsTag 

pattern  Plain       3    "Linear"  { 
puts $fileModel "pattern  Plain       3  Linear {"
	# Load    nodeTag    LoadValues 

    # SP    nodeTag    dofTag    DispValue 
    for { set j 0 } { $j <= $numY} {incr j 1} {
      #load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      #puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #load on 1 direction
      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #sp [expr $nodeYdirnStart+($numY+1)*$numZ+$j] 1 1.0
      #puts $fileModel "sp [expr $nodeYdirnStart+($numY+1)*$numZ+$j] 1 1.0"
    }
    for { set j 1 } { $j <= $numX} {incr j 1} {
      #load [expr 1+($numX+1)*$numZ+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      #puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000 +1.000000E+000   0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #load on 1 direction
      load [expr 1+($numX+1)*$numZ+$j]  +1.000000E+000 +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +1.000000E+000  +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
      #sp [expr 1+($numX+1)*$numZ+$j] 1 1.0 
      #puts $fileModel "sp [expr 1+($numX+1)*$numZ+$j] 1 1.0"
    }

    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
 
    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
    puts $fileModel "}"
} 

#pattern  Plain       2    "Linear"  { 
#puts $fileModel "pattern  Plain       2  Linear {"
#    # Load    nodeTag    LoadValues 
#    for { set j 0 } { $j <= $numY} {incr j 1} {
#      # load on 2 direction
#      #load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      #puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +0.000000E+000 +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#      #load on 1 direction
#      load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      puts $fileModel "load [expr $nodeYdirnStart+($numY+1)*$numZ+$j] +1.000000E+000 +0.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#    }
#    for { set j 1 } { $j <= $numX} {incr j 1} {
#      #load on 2 direction
#      #load [expr 1+($numX+1)*$numZ+$j]  +0.000000E+000  +1.000000E+000 0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      #puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +0.000000E+000 +1.000000E+000   0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#      #load on 1 direction
#      load [expr 1+($numX+1)*$numZ+$j]  +1.000000E+000 +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000
#      puts $fileModel "load [expr 1+($numX+1)*$numZ+$j] +1.000000E+000  +0.000000E+000  0.000000E+000  +0.000000E+000 +0.000000E+000 0.000000E+000"
#    }
#    # eleLoad    eleTags    beamUniform    Wy    Wz    <Wx> 
# 
#    # eleLoad    eleTags    beamPoint    Py    Pz    xL    <Px> 
#puts $fileModel "}"
#} 
close $fileModel 
set  r_1 [recorder  Node  -file  output_Node_DefoShape_Dsp.dis  -time -nodeRange 1  $endNodeTag -dof  1  2  3 disp]
#recorder  Element -file output_Ele_All.stress -time -eleRange 1 $endEleTag stress
#recorder  Element -file output_Ele_All.strain -time -eleRange 1 $endEleTag strain
set  r_2 [recorder Node -xml output_disp_11.xml -nodeRange 1  $endNodeTag -time -dof 1 2 3 4 5 6 disp]
#set  r_2 [recorder Node -file output_accel_11.txt -nodeRange 1  $endNodeTag  -time -dof 1 2 accel]
set  r_4 [recorder Element -xml output_stress_1_11.xml -time  -eleRange 1 $endEleTag stresses]
#set  r_5 [recorder Element -file output_stress_2_11.txt -time  -eleRange 1 $endEleTag material 2 forces]
#set  r_6 [recorder Element -file output_stress_3_11.txt -time  -eleRange 1 $endEleTag material 3 forces]
#set  r_7 [recorder Element -file output_stress_4_11.txt -time  -eleRange 1 $endEleTag material 4 forces]
set  r_8 [recorder Element -xml output_strain_1_11.xml -time  -eleRange 1 $endEleTag strains]
#set  r_9 [recorder Element -file output_strain_2_11.txt -time  -eleRange 1 $endEleTag material 2 deformations]
#set r_10 [recorder Element -file output_strain_3_11.txt -time  -eleRange 1 $endEleTag material 3 deformations]
#set r_11 [recorder Element -file output_strain_4_11.txt -time  -eleRange 1 $endEleTag material 4 deformations]
# ---------------------------------- 
# perform the analysis
# ---------------------------------- 
set IDctrlNode [expr 1+($numX+1)*$numZ];
set IDctrlDOF 1;
set Dincr 0.1;
# Constraint Handler Lagrange Transformation
constraints  Plain
# Convergence Test 
test  NormDispIncr  +1.000000E-002    250     0
#test EnergyIncr  1.0e-12    10         0
# Solution Algorithm 
algorithm  KrylovNewton 
# Integrator
integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
puts "integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr"
# DOF Numberer 
numberer  Plain 
# System of Equations 
system  SparseGeneral 
# Analysis Type 
analysis  Static 

set rcdcmd "recorder  Node  -file  PushoverCase_KeyNode_Dsp.out  -time -node [expr 1+($numX+1)*$numZ] -dof $IDctrlDOF  disp "
eval $rcdcmd;

#  perform Static Cyclic Displacements Analysis   /////////////////////////////////////////////////////////////

set iDmax1 [expr 0.08*$in/$mm];  # 2.032mm
set iDmax2 [expr 0.12*$in/$mm];  # 3.048mm
set iDmax3 [expr 0.3*$in/$mm];  #  7.62mm
set iDmax4 [expr 0.4*$in/$mm];  # 10.16mm
set iDmax5 [expr 0.6*$in/$mm];  # 15.24mm
set iDmax6 [expr 1.1*$in/$mm];  # 27.94mm

# perform the analysis    /////////////////////////////////////////////////////////////
#                0.08   -0.08  0.12   -0.12   0.3   -0.3   0.4   -0.4   0.6   -0.6   1.1   -1.1
#set numSteps     {   508   1016  1270    1524  2667   3810  4445   5080  6350   7620 10795  13970 }
set numIters     {     0    200   400     400   400    400   400    400   400    400   400    400 }
#set increments   { 0.004 -0.004 0.004  -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 0.004 -0.004 }
set numSteps     { 20   40   50   60  106  152 178  204 254  304 432  560 }
set increments   { 0.1 -0.1 0.1  -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 }
	
for { set i 0 } { $i<12 } { incr i 1 } {
    set numStep [lindex $numSteps $i]
    set numIter [lindex $numIters $i]
    set Dincr [lindex $increments $i]
    integrator DisplacementControl $IDctrlNode $IDctrlDOF $Dincr
    if {$numIter == 0} {
	set numIter 100
    }
    test NormDispIncr 1.0e-3 $numIter 5
    #analyze $numStep

    set j 0;
    while { $j < $numStep } {
      set ok [analyze 1]
      # ----------------------------------------------if convergence failure-------------------------
      # if analysis fails, we try some other stuff
      # performance is slower inside this loop  global maxNumIterStatic;      # max no. of iterations performed before "failure to converge" is ret'd
      if {$ok != 0} {
        puts "Trying Newton with Initial Tangent .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm Newton -initial
        set ok [analyze 1]
        test NormDispIncr 1e-3 [expr $numIter*2] 5
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying Broyden .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm Broyden 8
        set ok [analyze 1 ]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        puts "Trying NewtonWithLineSearch .."
        test NormDispIncr 1e-3 $numIter 0
        algorithm NewtonLineSearch 0.8 
        set ok [analyze 1]
        algorithm KrylovNewton
      }
      if {$ok != 0} {
        set putout "PROBLEM: $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF]";
        puts $putout
        test NormDispIncr 1e-3 $numIter 5
        set ok [analyze 1]
      }; # end if
      incr j 1;
      #record;
      #puts $end
    }

    puts $i

}
#remove recorders;

