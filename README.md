# AccelerateFFT
Simple real FFT using Mac OS X Accelerate

create a command line tool, for eexample called "AccelerateTool". click on the blue Xcode porject icon (the name of project, eg "AccelerateTool") and then add the Accelerate framework:

Go into Build Phases 
Link Binary with Libraries
choose Acclerate

the example wraps the code form Apple's documentation
https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-SW1

to do the FFT and return magnitude



