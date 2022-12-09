# EDMA_FDS_analysis_2D

Matlab code for computing Facial Dismorphology Score (FDS) via iterative bootstrapping to assess the statistical significance of the patterns revealed by EDMA.

This function generates random groups (with an increasing number of patients at each iteration) and computes the percentage of significant inter-landmark differences against a random control grup.

As an output, a FDS histogram for each iteration is generated and saved as a figure.

%%%%%%%%%%%%%%%%%%%%
% Input:
% res_path: Result path where the histograms will be stored.
% nameOutput: name for the stored histogram.
% id_cmp: 1 for 'Down', 2 for 'Morquio', 3 for 'NFM' (neurofibromatosis), and 4 for 'Noonan'
% age_range: [min max] vector defining the range of ages.
% num_samples: number of samples in the random groups.
% sample_step: number of patient added in each iteration of the bootstrap experiment
%
%
% Example:
% EDMA_FDS_analysis_2D('results\', 'Morquio FDS Histogram', 2, [0 14], 9, 3)
%%%%%%%%%%%%%%%%%%%%

This code is available as Supplementary Material to the article "Population-specific facial dysmorphologies and diagnosis accuracy of genetic and rare diseases in an admixed Colombian population" authored by: Luis Miguel Echeverry, Estephania Candelo, Eidith Gómez, Paula Solís, Diana Ramírez, Diana Ortiz, Alejandro González, Xavier Sevillano, Juan Carlos Cuéllar, Harry Pachajoa, and Neus Martínez-Abadías
