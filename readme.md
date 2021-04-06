# ECE-GY 6023.  Introduction to Wireless Communications

This repository provides instructional material for the
graduate wireless communications, ECE-GY 6023, at New York University
taught by [Sundeep Rangan](http://wireless.engineering.nyu.edu/sundeep-rangan/):

Anyone is free to use and copy this material (at their own risk!).
But, please cite the material if you use the material in your own class.

## Pre-requisites

The course assumes you are familiar with digital communications at the graduate level.  There are many resources for digital communications, including some lecture notes I created for the [NYU class](https://github.com/sdrangan/digitalcomm).

Additionally, some lecture notes (and problems to be added later) assume you have access to MATLAB along with the communications, phased array and antenna toolboxes.

## Feedback

Any feedback is welcome.  If you find errors, have ideas for improvements,
or want to voice any other thoughts, [create an issue](https://help.github.com/articles/creating-an-issue/)
and we will try to get to it.
Even better, fork the repository, make the changes yourself and
[create a pull request](https://help.github.com/articles/about-pull-requests/)
and we will try to merge it in.  See the [excellent instructions](https://github.com/ishjain/learnGithub/blob/master/updateMLrepo.md)
from the former TA Ish Jain.


## Lecture Sequence
The tentative plan for the lectures are below.  Right now, only a few lectures
have full material.  We will be hoping to add to this material over the course
of the semester.  Other topics may be added at the end depending on time.

* Course Introduction
    * Lecture: [[PDF]](./lectures/CourseAdmin.pdf) [[PPT]](./lectures/CourseAdmin.pptx) 
    * Lecture video:  [[YouTube]](https://youtu.be/DZLp12GCHow)
* Unit 1.  Basics of Antennas and Free-space Propagation 
    * Lecture: [[PDF]](./lectures/Unit01_Antennas.pdf) [[PPT]](./lectures/Unit01_Antennas.pptx) 
    * [Lecture videos](./unit01_antennas/readme.md) and in-class exercises
    * Demo: Calculating and displaying antenna patterns [[PDF]](./unit01_antennas/demo_antennas.pdf) [[Matlab]](./unit01_antennas/demo_antennas.m)
    * Demo: 3GPP 5G antenna model [[PDF]](./unit01_antennas/demo_3gpp_antenna.pdf) [[Matlab live]](./unit01_antennas/demo_3gpp_antenna.m)
    * Problems:  [[PDF]](./unit01_antennas/prob/prob_antennas.pdf) [[Latex]](./unit01_antennas/prob/prob_antennas.tex)
    * Lab:  Simulating a 28 GHz antenna for a UAV [[PDF]](./unit01_antennas/lab_uav_antenna_partial.pdf) [[Matlab]](./unit01_antennas/lab_uav_antenna_partial.m)
* Unit 2.  Non-LOS Propagation and Link-Budget Analysis 
    * Lecture: [[PDF]](./lectures/Unit02_Propagation.pdf) [[PPT]](./lectures/Unit02_Propagation.pptx) 
    * [Lecture videos](./unit02_propagation/readme.md) and in-class exercises
    * Demo: Simulating AWGN Noise [[PDF]](./unit02_propagation/demo_awgn.pdf) [[Matlab Live]](./unit02_propagation/demo_awgn.mlx)
    * Demo: Propagation and rate modeling [[PDF]](./unit02_propagation/demo_path_loss_model.pdf) [[Matlab]](./unit02_propagation/demo_path_loss_model.m)
    * Problems:  [[PDF]](./unit02_propagation/prob/prob_propagation.pdf) [[Latex]](./unit02_propagation/prob/prob_propagation.tex)
    * Lab:  Propagation modeling from ray tracing data [[PDF]](./unit02_propagation/lab_prop_modeling_partial.pdf) [[Matlab]](./unit02_propagation/lab_prop_modeling_partial.m)
* Unit 3.  Multipath Fading
    * Lecture: [[PDF]](./lectures/Unit03_Fading.pdf) [[PPT]](./lectures/Unit03_Fading.pptx) 
    * [Lecture videos](./unit03_fading/readme.md) and in-class exercises
    * Demo: Simulating fading [[PDF]](./unit03_fading/demo_fading.pdf) [[Matlab Live]](./fading/unit03_demo_fading.mlx)
    * Problems:  [[PDF]](./unit03_fading/prob/prob_fading.pdf) [[Latex]](./unit03_fading/prob/prob_fading.tex)
    * Lab:  5G channel sounding with Doppler [[PDF]](./unit03_fading/partial/lab_chan_sounder.pdf) [[Matlab]](./unit03_fading/partial/lab_chan_sounder.m)
* Unit 4.  Capacity and Coding on Fading Channels
    * Lecture: [[PDF]](./lectures/Unit04_Coding.pdf) [[PPT]](./lectures/Unit04_Fading.pptx) 
    * [Lecture videos](./unit04_coding/readme.md) and in-class exercises
    * Demo: Uncoded BER on fading channel [[PDF]](./unit04_coding/demo_uncoded.pdf) [[Matlab Live]](./unit04_coding/demo_uncoded.mlx)
    * Demo: Convolutional coding on a fading channel[[PDF]](./unit04_coding/demo_conv.pdf) [[Matlab Live]](./unit04_coding/demo_conv.mlx)
    * Lab:  5G NR Downlink Throughput with Fading and LDPC coding [[PDF]](./unit04_coding/lab_partial/labPdsch.pdf) [[Matlab Live]](./unit04_coding/lab_partial/labPdsch.mlx)
    * Problems:  [[PDF]](./unit04_coding/prob/prob_coding.pdf) [[Latex]](./unit04_coding/prob/prob_coding.tex)
* Unit 5.  Adaptive Modulation and Coding
    * Lecture: [[PDF]](./lectures/Unit05_AMC.pdf) [[PPT]](./lectures/Unit05_AMC.pptx) 
    * [Lecture videos](./unit05_amc/readme.md) and in-class exercises
    * Demo: 802.11 MCS selection  [[PDF]](./unit05_amc/demo_mcs.pdf) [[Matlab Live]](./unit05_amc/demo_mcs.mlx)
    * Demo: Channel Tracking with 5G NR CSI-RS [[PDF]](./unit05_amc/demo_csirs.pdf) [[Matlab Live]](./unit05_amc/demo_csirs.mlx)    
    * Lab:  5G NR DL Throughput with Multi-Process HARQ [[PDF]](./unit05_amc/lab_partial/labHarq.pdf) [[Matlab Live]](./unit05_amc/lab_partial/labHarq.mlx) 
    * Problems:  [[PDF]](./unit05_amc/prob/prob_amc.pdf) [[Latex]](./unit05_amc/prob/prob_amc   .tex)
* Unit 6.  Diversity
* Unit 7.  OFDM Channel Estimation and Equalization
    * Lecture:  [[PDF]](./lectures/Unit07_ChanEst.pdf) [[Powerpoint]](../lectures/Unit07_ChanEst.pdf) 
    * [Lecture videos](./unit07_chanest/readme.md) and in-class exercises
    * Demo:  5G NR DM-RS configuration  [[PDF]](./unit07_chanest/demoDMRSConfig.pdf)  [[Matlab Live]](./unit07_chanest/demoDMRSConfig.mlx) 
    * Demo:  Kernel regression channel estimation [[PDF]](./unit07_chanest/demoKernelEst.pdf)  [[Matlab Live]](./unit07_chanest/demoKernelEst.mlx) 
    * Lab:  5G NR DL Throughput with Channel Estimation [[PDF]](./unit07_chanest/lab_partial/labChanEst.pdf) [[Matlab Live]](./unit07_chanest/lab_partial/labChanEst.mlx) 
* Unit 8.  Multiple Antennas and Beamforming
    * Lecture: [[PDF]](./lectures/Unit06_Beamforming.pdf) [[PPT]](./lectures/Unit06_Beamforming.pptx) 
    * Demo: Visualizing and simualting arrays [[PDF]](./beamforming/demo_bf.pdf) [[Matlab]](./beamforming/demo_bf.m)
    * In-class exercises: [[PDF]](./beamforming/bf_inclass_partial.pdf) [[Matlab]](./beamforming/bf_inclass_partial.m)
    * Lab:  5G NR downlink simulation with beamforming [[PDF]](./beamforming/partial/lab_pdsch_bf.pdf) [[Matlab]](./beamforming/partial/lab_pdsch_bf.m)
* Unit 9.  Beam Tracking and Directional Estimation
* Unit 10.  Introduction to MIMO 


