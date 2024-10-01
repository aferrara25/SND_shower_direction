# Function used to combine plots related to the same energy from different run
cd ./output
hadd -f 3W300Gev.root TB_outputRun_100642.root TB_outputRun_100643.root TB_outputRun_100639.root 
hadd -f 3W240Gev.root TB_outputRun_100646.root TB_outputRun_100647.root TB_outputRun_100648.root
hadd -f 3W180Gev.root TB_outputRun_100671.root TB_outputRun_100672.root TB_outputRun_100635.root TB_outputRun_100636.root
hadd -f 3W140Gev.root TB_outputRun_100633.root TB_outputRun_100673.root TB_outputRun_100674.root
hadd -f 3W100Gev.root TB_outputRun_100630.root TB_outputRun_100632.root 