%% DCE MRI Analysis Code to generate LRRM and Tofts analysis and their parametric maps
%This code works only acquired from Bruker preclinical scanner (invivo or invitro).
%This code only works if acquired anatomical scans with RAREVTR sequence for T1 maps
  As it also generates T1 parametric maps.
%This code only works if acquired DCE-MRI scans with FLASH sequence
%In our sequence, we have use two different scans for Pre and Post FLASH
%For Pre SCE-MRI scans, we have used 35 time frames with 10 s temporal resolution.
%For Post DCE-MRI scans, we have used 135 time frames with 10 s temporal resolution.
%We have used gadavist as a MRI contrast agent, if you use another contrast agent please change the R1 of that contrast agent, you can change in the line 181 in dce_mri.m
%If changing the time frames or the temporal resolution of FLASH sequence, need to change the 	       lines in the dce_mri.m from lines 45-48, 326 & 338 according to your time frames (as it is hard coded!!!) and line 341 according to your temporal resolution.
