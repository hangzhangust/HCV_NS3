%%
% *Epistatic interactions promote resistance against 
% direct-acting antivirals for HCV NS3 protein*

% Code for developing the model and re-generating the figures in the paper

%% Setting up paths (of functions and data files required) and necessary parameters

clear;
close all;
clc;

addpath(genpath("data"));
addpath functions
addpath 3rd_party_code
addpath figures_scripts

%% Data preparation and model training

% HCV NS3 1a sequences were downloaded from GLUE data base
% (http://hcv-glue.cvr.gla.ac.uk; accessed Sep. 2, 2021)
% load sequences, patient id, outliers and sequence id without patient information
load('data_NS3.mat')


% remove outliers  and sequence without information;
remove_id = unique([No_info;outliers]);
header(remove_id)=[];
sequences(remove_id)=[];
patient_id(remove_id)=[];


% re-weight sequences according to number of sequences from each patient
weight_seq = get_weight_seq(patient_id);

% feed the sequences and weight_seq to the MPF-BML-GUI to train the model
% Code for running MPF-BML is freely available at <https://github.com/ahmedaq/MPF-BML-GUI>. 

%% Fig. 1 Validation of the inferred NS3 fitness landscape.

% plot the energy vs fitness of the inferred model and compare with the
% conservation-only model
clear;
clc;
run figure1.m

%% Fig. 2 SC-DRMs are enriched in NS3 drugs.
% (a) list of NS3 drugs and the associated DRMs. check data/Drugs_mutants.mat
% (b) statistical significance of the number of SC-DRMs associated with each drug.
clear;
clc;
run figure2.m


%% Fig. 3 Identification of SC-DRMs and their significance.

% (a) network of interactions between top 10/100/300 ranked mutations 
% (b) pairs of interacting residues involving SC-DRMs that are in contact 
% based on the crystal structure of the NS3 protein. check data/4b6e.pse file
% (c) inferred NS3 sectors and their association with SC-DRMs. Code for 
% running RocaSec to get the sectors is freely available at 
% <https://github.com/ahmedaq/RocaSec>. 
clear;
clc;
run figure3.m

%% Fig. 4 Histogram of the change in energy observed by all single mutations X in 
% the H77 strain carrying
% (a) the D168E mutant and (b}) the Q80K mutant. 

clear;
clc;
run figure4.m
%% Fig. 5 SC-DRMs appear to directly impact binding of NS3 drugs.
% (a) binding residues of drugs shown on the crystal structure of the NS3 protein-drug complexes
% check data/*.pse file (* can be danoprevir/grazoprevir/telaprevir/vaniprevir)
% (b-c) statistical significance of the number of (b) drug-specific DRMs/SC-DRMs and 
% (c) all DRMs/SC-DRMs in binding residues of each of the four considered drugs.
clear;
clc;
run figure5.m

%% Fig. 6 Escape time of residues involved in NS3 DRMs.
% (a) comparison between escape time of residues involved in SC-DRMs 
% and the remaining residues involved in DRMs.
% (b) individual escape time of residues involved in DRMs of the NS3 protein.
clear;
clc;
run figure6.m

%% Fig. 7 Correlation between the drug efficacy and the number of SC-DRMs associated with each drug.
% (a) NS3-specific DAAs and (b) multi-protein DAAs.
clear;
clc;
run figure7.m
%% Supp. Fig. 1 Robustness of the correlation observed between model energies
% and experimental fitness values. Results are shown for the maximum-entropy model 
% that considers epistatic interaction and 
% for the conservation-only model that ignores epistasis.
clear;
clc;
run suppfigure1.m


%% Supp. Fig. 2 Top ranked pairs of mutations based on the strength of their couplings.
clear;
clc;
run suppfigure2.m

%% Supp. Fig. 3 Robustness of the enrichment of SC-DRMs in each drug to the number of top-coupled pairs of mutations used to define SC-DRMs.
clear;
clc;
run suppfigure3.m

%% Supp. Fig. 4 Statistical significance of the number of non-SC-DRMs associated with each drug.
clear;
clc;
run suppfigure4.m
%% Supp. Fig. 5 Comparison of statistical properties and model predictions 
% based on complete data with those based on a subset of drug-naïve patients.
% (a) Correlation of single mutant probabilities (left panel) and double mutant
% probabilities (right panel) between sequences from all patients (7370 sequences) 
% and the subset of drug-naïve patients (5877 sequences).
% (b) Correlation between model predicted energies and experimental 
% fitness measurements compiled from different studies.

clear;
clc;
run suppfigure5.m
%% Supp. Fig. 6 Correlation between the model predicted energies and 36 experimental 
% fitness measurements that are associated with DRMs.

clear;
clc;
run suppfigure6.m
%% Supp. Fig. 7 Correlation between the model predicted energies and 36 experimental 
% fitness measurements that are associated with DRMs.

clear;
clc;
run suppfigure7.m

%% Supp. Fig. 8 Statistical validation of the inferred HCV NS3
clear;
clc;
run suppfigure8.m
