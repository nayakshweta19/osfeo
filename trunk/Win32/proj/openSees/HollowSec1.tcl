# HollowSec1.tcl 

######################## 
# Analysis-Sequence  1 #
######################## 

# Start of model generation 
# ========================= 

# Create ModelBuilder 
# ------------------- 
model  BasicBuilder  -ndm  2  -ndf  3 

# Define geometry 
# --------------- 
source  NodeCoord.tcl 

# Define Single Point Constraints 
# ------------------------------- 
source  SPConstraint.tcl 

# Define nodal masses 
# ------------------- 
source  NodeMass.tcl 

# Define Multi Point Constraints 
# ------------------------------ 
source  MPConstraint.tcl 

# Define material(s) 
# ------------------ 
source  Materials.tcl 

# Define section(s) 
# ----------------- 
source  Sections.tcl 

# Define geometric transformation(s) 
# ---------------------------------- 
source  GeoTran.tcl 

# Define element(s) 
# ----------------- 
source  Elements.tcl 

# Define time series 
# ------------------ 
source  TimeSeries.tcl 

# Start of anaysis generation 
# =========================== 

# Get Initial Stiffness 
# --------------------- 
initialize 

# Analysis: StaticDefaultCase 
# +++++++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
source  LoadPattern_1.tcl 

# Define recorder(s) 
# -------------------- 
source  Recorder_1.tcl 

# Define analysis options 
# ----------------------- 
source  AnalysisOptn_3.tcl 

analyze     1 

# Reset for next analysis case 
# ---------------------------- 
loadConst -time 0.0 
remove recorders 
wipeAnalysis 

# Analysis: EigenDefaultCase 
# ++++++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
source  LoadPattern_2.tcl 

# Define recorder(s) 
# -------------------- 
source  Recorder_2.tcl 

# Define analysis options 
# ----------------------- 
source  AnalysisOptn_2.tcl 

set eigFID [open EigenDefaultCase_Node_EigenVector_EigenVal.out w] 
puts $eigFID [eigen  generalized  fullGenLapack     4] 
close $eigFID 
analyze  1  0.0001 

# Reset for next analysis case 
# ---------------------------- 
setTime 0.0 
remove recorders 
wipeAnalysis 

# Analysis: PushoverCase 
# ++++++++++++++++++++++ 

# Define load pattern 
# ------------------- 
source  LoadPattern_3.tcl 

# Define recorder(s) 
# -------------------- 
source  Recorder_3.tcl 

# Define analysis options 
# ----------------------- 
source  AnalysisOptn_3.tcl 

analyze    30 

# Clean up 
# -------- 
wipe 
exit 

