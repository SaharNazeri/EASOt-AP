# EASOt-AP
EASOt-AP: An open-source MATLAB code package to Estimate the Earthquake Source Properties: Seismic Moment, Rupture Radius and Stress-drop from Time Evaluation of P-wave Displacement Amplitude

Corresponding Author (sahar.nazeri@unina.it)
Department of Physics E. Pancini, University of Naples Federico II, Naples, Italy

Main Reference: Zollo, A., Nazeri, S., Colombelli, S., 2021. Earthquake Seismic Moment, Rupture Radius, and Stress Drop From P-Wave Displacement Amplitude Versus Time Curves. IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-11, 2022, Art no. 5912211, doi: 10.1109/TGRS.2021.3119909

We developed a Matlab package, to rapidly obtain the apparent EArthquake SOurce time function using the Average of near-source P-wave content of the observations (EASOt-AP) and estimate the earthquake source parameters, such as the seismic moment, the rupture radius, and the average static stress drop. The algorithm implemented in this package is based on a rapid and straightforward methodology that is recently developed by Zollo et al. (2021). To this purpose, EASOt-AP retrieves and models the time evolution of the average P-wave displacement in the logarithm scale named as LPDT curve. In this paper, the performance of the EASOt-AP has been shown by evaluating the strong motion data recorded in the near-source range of small magnitude events (2≤Mw≤3) that occurred in the Irpinia region, southern Italy. This Tool is designed for an easy use, and all steps of data processing are performed completely within the numerical environment of MATLAB.  
Description of the software, Algorithm, and Implementation
The EASOt-AP package is designed using Matlab (R2021a) supported on both operating Systems of Windows and Linux/Ubuntu, thus it is recommended to run it using the same version or later releases of Matlab. The design philosophy of the package divides it into 3 main parts of input (folder “INPUT”), analysis box (folder “CODES”) and output (folder “OUTPUT”). 
Following section enumerate all input parameters in more details. 

Input
This package requires a SAC formatted strong motion data with the following needed fields in the header: (1) earthquake location (longitude, latitude, depth) and origin time, (2) station coordinates (longitude, latitude, elevation), (3) P-wave arrival time, stored as “A” field, (4) source receiver distance in kilometre stored as “DIST” field (automatically calculated by SAC after setting the earthquake and station locations), (5) magnitude of the event (M_w is preferred) reported by agencies as a piece of pre-information (stored as “MAG” field). 
It is expected that a quality check of the waveforms is done before running the package, for instance, whenever manually or automatically picking the P-phase. However, the analysis module is provided with checking the quality of data with signal-to-noise ratio (SNR) larger than a threshold value assigned by user as an input parameter. Also, it is important to use of the waveforms with at least 2 seconds pre-event signal, otherwise the waveforms are excluded automatically from the calculation box.
In addition to the strong motion data, the main input file used by EASOt-AP package is a text file named as “Input.txt” saved in the folder of “INPUT”. To increase the flexibility of the data analysis step, this file is provided with a variety of optional parameters/arguments that should be modified for any new analyses or datasets. Description of the set-up of these inputs and their default values are presented as following:
(1) “Region”: name of the studied area, same as the name of the main directory of data i.e., (Package/INPUT/”Region”). Inside this folder, SAC files of all events are presented in separate subfolders for each event that can be specified with any kind of name, for example date and time of the occurrence of the event. After running the package, a directory with the same name of the region will be created automatically inside the folder of the “OUTPUT” as (Package/OUTPUT/”Region”).
(2)  “wave-type”: type of the raw ground motion data as “A” (Acceleration) or “V” (Velocity). (Default: A).
(3) “dirparam”: specific part of the name of the SAC-formatted waveforms to make directory. (Default: *sac*).
(4) “CF, Unit”: CF is a conversion factor to change the conventional unit of amplitude to the SI unit i.e., m/s2 or m/s for acceleration or velocity respectively. “Unit” can be set as Yes or No depends on existence of a unique CF value for whole data or not.  In the case of setting “Unit” as “No”, different conversion factors corresponding to each waveform will be read from the header of the SAC formatted data stored as the “scale” field. (Default: 1, No). 
(5) “minSta”: minimum number of stations required to make average and build up the LPDT curve at each time sample. (Default: 5)
(6) “Rmax”: maximum hypocentral distance in km to select data within, for events with different range of Magnitude: M?1, 1<M?3, 3<M?5, and M>5. (Default: 20, 40, 60, 100).
(7) “SNR”, the threshold value for SNR to exclude the bad quality signals. (Default: 10).
(8) “nPol, cor1, cor2”: number of poles (nPol) and two corner frequencies (cor1, cor2 in Hz) of the band-pass Butterworth filter applied to the displacement. (Default: 2, 1, 30).
(9) “Vfilter, nPolV, cor1V, cor2V”: number of poles (nPolV) and two corner frequencies (cor1V, cor2V in Hz) of the band-pass Butterworth filter applied to the velocity signals in case that “Vfilter” is set as “Yes”. If “Vfilter” is assigned as “No”, no filter is applied to the velocity waveforms. In the other words, this step is optional based on the “Vfilter” value. (Default: No, 2, 0.075, 30).
(10) “WinFix, WinMax, S_P_Coef”: the P-waves are windowed after the onset using a fix (WinMax) or variable values for each station, depending on what is set as the “WinFix”, i.e., “Yes” or “No” respectively. In case of selecting the variable time window, the S-mines-P coefficient (S_P_Coef) is then used to calculate the length of the window from the P-wave to the S-wave. (Default: No, 1, 0.09).
(11) “WinFix_fit, WinMax_fit”: selecting a fix time window (WinFix_fit set as “Yes”) to fit to the LPDT curve with length of “WinMax_fit”. If “WinFix_fit” is assigned as “No”, the entire of the LPDT curve is used to find the best nonlinear fit to the curve. (Default: No, 1). 
(12) “rho, Vp, Vs, SD”: properties of the medium i.e., density of the area (rho in kg/m3), P-wave velocity (Vp in m/s), S-wave velocity (Vs in m/s), and trial value of Stress Drop (SD in MPa). (Default: 2700, 5500, 3000, 1)
(13) “RadP”: P-wave radiation pattern coefficient. (Default: 0.52). 
(14) “QFilter, Qp”: implementing (on) or ignoring (off) the anelastic attenuation analysis. In case of implementing this correction, the constant Qp-value for the studied area is assumed to remove the anelastic attenuation from the final solution of STF. (Default: on, 150).
(15) “FitFun”, selecting the fit function to model the LPDT curve, among two options of “exp” or “HB”, (see section C. Modelling the LPDT curve and Table 1 for more detail). (Default: “HB”).
(16) “SD1, SD2, SDN”: range of the static stress drop in MPa from SD1 to SD2 which this interval is divided to SDN values to find the best stress drop value implementing the “Approach 1”.

How to run the package:
To run the package, follow the below structure:
1.	create a main folder with the name of the package as “EASOt-AP”.
2.	Create two subfolders as “.../EASOt-AP/CODES”, and “.../EASOt-AP/INPUT”.
3.	Download all materials and copy pate to the relevant folders.
4.	Open the main file in the directory of CODES named as “SaharLPXT_Main.m”. 
5.	Change the current path of Matlab to the path of the package on the pc as: “…/../EASOt-AP”.
6.	Run the code.
7.	Check the result if the folder of “../EASOt-AP/OUTPUT”.
