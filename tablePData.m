% Script that organizes all of the data for viewing and QC'ing.
clear;

T=struct('ID',[],'file',[],'Mc',[],'Tf',[]);

T(1).ID='BnE_c4';
T(1).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/BnE_c4.mat';
T(1).Mc=1.5;
T(1).dM=0.2;
T(1).Tf=[0 Inf];

T(2).ID='BnE_c5';
T(2).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/BnE_c5.mat';
T(2).Mc=1.5;
T(2).dM=0.2;
T(2).Tf=[0 Inf];

T(3).ID='ESB07';
T(3).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/ESB07.mat';
T(3).Mc=0.0;
T(3).dM=0.2;
T(3).Tf=[0 Inf];

T(4).ID='ESB08';
T(4).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/ESB08.mat';
T(4).Mc=0.2;
T(4).dM=0.2;
T(4).Tf=[0 Inf];

T(5).ID='ESB09';
T(5).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/ESB09.mat';
T(5).Mc=0.0;
T(5).dM=0.2;
T(5).Tf=[0 Inf];

T(6).ID='ESB10';
T(6).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/ESB10.mat';
T(6).Mc=0.6;
T(6).dM=0.2;
T(6).Tf=[0 Inf];

T(7).ID='Pena';
T(7).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Pena.mat';
T(7).Mc=0.8;
T(7).dM=0.2;
T(7).Tf=[0 Inf];

T(8).ID='PH';
T(8).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/PH.mat';
T(8).Mc=-5.0;
T(8).dM=0.2;
T(8).Tf=[0 Inf];

T(9).ID='PNR1z';
T(9).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/PNR1z.mat';
T(9).Mc=-1.0;
T(9).dM=0.1;
T(9).Tf=[0 Inf];

T(10).ID='PNR2';
T(10).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/PNR2.mat';
T(10).Mc=-1.0;
T(10).dM=0.1;
T(10).Tf=[0 Inf];

T(11).ID='AR-c1';
T(11).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c1.mat';
T(11).Mc=0.2;
T(11).dM=0.1;
T(11).Tf=[0 Inf];

T(12).ID='AR-c6';
T(12).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c6.mat';
T(12).Mc=0.2;
T(12).dM=0.1;
T(12).Tf=[0 Inf];

T(13).ID='AR-c7';
T(13).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c7.mat';
T(13).Mc=0.0;
T(13).dM=0.1;
T(13).Tf=[0 Inf];

T(14).ID='AR-c8';
T(14).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c8.mat';
T(14).Mc=0.2;
T(14).dM=0.1;
T(14).Tf=[0 Inf];

T(15).ID='AR-c11';
T(15).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c11.mat';
T(15).Mc=0.5;
T(15).dM=0.1;
T(15).Tf=[0 Inf];

T(16).ID='AR-c15';
T(16).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c15.mat';
T(16).Mc=0.2;
T(16).dM=0.1;
T(16).Tf=[0 Inf];

T(17).ID='AR-c16';
T(17).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Ark_c16.mat';
T(17).Mc=0.0;
T(17).dM=0.1;
T(17).Tf=[0 Inf];

T(18).ID='SSFS-1993';
T(18).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/SSFS1993.mat';
T(18).Mc=0.1;
T(18).dM=0.5;
T(18).Tf=[0 Inf];

T(19).ID='SSFS-2000';
T(19).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/SSFS2000.mat';
T(19).Mc=0.0;
T(19).dM=0.5;
T(19).Tf=[0 Inf];

T(20).ID='SSFS-2003';
T(20).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/SSFS2003.mat';
T(20).Mc=-0.2;
T(20).dM=0.4;
T(20).Tf=[0 Inf];

T(21).ID='SSFS-2004';
T(21).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/SSFS2004.mat';
T(21).Mc=-0.5;
T(21).dM=0.4;
T(21).Tf=[0 Inf];

T(22).ID='SSFS-2005';
T(22).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/SSFS2005.mat';
T(22).Mc=-0.5;
T(22).dM=0.4;
T(22).Tf=[0 Inf];

T(23).ID='CB-HAB1';
T(23).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/CB_HAB1.mat';
T(23).Mc=-0.5;
T(23).dM=0.2;
T(23).Tf=[0 Inf];

T(24).ID='CB-HAB4';
T(24).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/CB_HAB4.mat';
T(24).Mc=1.3;
T(24).dM=0.4;
T(24).Tf=[0 Inf];

T(25).ID='Gross Schoenbeck';
T(25).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/GS.mat';
T(25).Mc=-1.3;
T(25).dM=0.2;
T(25).Tf=[0 Inf];

T(26).ID='StGallen';
T(26).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/StGallen.mat';
T(26).Mc=-0.5;
T(26).dM=0.1;
T(26).Tf=[0 Inf];

T(27).ID='Basel';
T(27).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Basel.mat';
T(27).Mc=0.7;
T(27).dM=0.3;

T(28).ID='Helsinki';
T(28).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Helsinki.mat';
T(28).Mc=0.2;
T(28).dM=0.1;
T(28).Tf=[0 Inf];

T(29).ID='Paralana';
T(29).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Paralana.mat';
T(29).Mc=0.3;
T(29).dM=0.3;
T(29).Tf=[0 Inf];

T(30).ID='KTB';
T(30).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/KTB.mat';
T(30).Mc=-0.7;
T(30).dM=0.1;
T(30).Tf=[0 Inf];

T(31).ID='Pohang';
T(31).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Pohang.mat';
T(31).Mc=0.5;
T(31).dM=0.2;
T(31).Tf=[0 Inf];

T(32).ID='BerlinGF';
T(32).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/BerlinGF.mat';
T(32).Mc=0.9;
T(32).dM=0.2;
T(32).Tf=[0 Inf];

T(33).ID='Youngstown';
T(33).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Youngstown.mat';
T(33).Mc=0.8;
T(33).dM=0.2;
T(33).Tf=[0 Inf];

T(34).ID='ParaV';
T(34).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/ParaV.mat';
T(34).Mc=0.5;
T(34).dM=0.3;
T(34).Tf=[0 Inf];

T(35).ID='GGB';
T(35).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/GGB.mat';
T(35).Mc=0.5;
T(35).dM=0.1;
T(35).Tf=[0 Inf];

T(36).ID='Castor';
T(36).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Castor.mat';
T(36).Mc=0.6;
T(36).dM=0.2;
T(36).Tf=[0 Inf];

T(37).ID='AzleReno';
T(37).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/AzleReno.mat';
T(37).Mc=0.6;
T(37).dM=0.15;
T(37).Tf=[0 Inf];

T(38).ID='Irving';
T(38).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Irving.mat';
T(38).Mc=2;
T(38).dM=0.15;
T(38).Tf=[0 Inf];

T(39).ID='OK';
T(39).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/OK.mat';
T(39).Mc=3.0;
T(39).dM=0.25;
T(39).Tf=[0 Inf];

T(40).ID='Rongchang';
T(40).file='/Users/rjs10/Desktop/Trailing/IS-bath/data/processed/Rongchang.mat';
T(40).Mc=2.0;
T(40).dM=0.2;
T(40).Tf=[0 Inf];
