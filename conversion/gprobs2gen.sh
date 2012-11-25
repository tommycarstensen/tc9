more +2 out_BEAGLE/BeagleOutput.X.bgl.ProduceBeagleInput.X.bgl.gprobs | awk -F : '{print $1, $2}' | awk '{print $1":"$2, $1":"$2, $0}'  | cut -d " " -f 3 --complement > out_BEAGLE/BeagleOutput.X.gen
