# Sonar3D
Code to process raw RW2 files produced by a "3D Profiling Sonar 2001" (Marine Electronics - http://www.marine-electronics.co.uk/3D_2001.html)

The UtrechtCoastalGroup owns a 3D Profiling Sonar 2001, manufactured by Marine Electronics. Although Marine Electronics provides code to process its RW2 files, we thought we could do a better job in detecting the bed. We developed this (Matlab) code within the EU-funded Hydralab IV Bardex II project. Results of the subsequent analyses can be found in Ruessink et al. (2015), https://doi.org/10.3390/jmse3041568. Later on, we were commissioned by Rijkswaterstaat to extend and apply the code to data collected as part of the Kustgenese 2 project. As part of this work, we made the resulting code OpenSource through the present GitHub. The code is now free for everyone to use.

This repository contains code in which bedlevel data are low-passed filtered with a quadratic loess filter. This code was written by Nathaniel Plant and has been documented in a paper in Marine Geology, please see https://doi.org/10.1016/S0025-3227(02)00497-8.

If you use the code in this repository, we would kindly ask you to refer to both Ruessink et al. (2015) - as the initial source of the code - and to Plant et al. (2002) - as the source for the loess code.

Please start with Process_sonar_example.m and have a look at the Word document (in Dutch).


