# Notebook for Amorim, Di et al. "Evolutionary Consequences of Domestication on the Selective Effects of New Amino Acid Changing Mutations in Canids"
# Author: Carlos Eduardo G. Amorim, guerraamorim[at]gmail
# Date: July 31, 2025

# you will need to have dadi-cli and varDFE set up
# https://dadi-cli.readthedocs.io/en/latest/
# https://github.com/meixilin/varDFE/

HERE=/path/to/your/main/directory
SFSDIR=/path/to/your/SFS/directory
VCFDIR=/path/to/your/VCF/directory # VCF files can be downloaded from here: https://datadryad.org/share/Xu8h6-72I-JWM1MFyKS2dMUwvIYYSiWV1SsgFHShJ3E

OUT=$HERE/$RUN

##########################################################################################
########################################## Compute SFSs ##################################
##########################################################################################

# Create SFSs with dadi-cli.
# Below are the commands used to generate the SFSs for the whole exome, considering the four high coverage (HC) population samples, and the projection values shown below.
# To generate the SFS for gene subsets, different population samples, or different projection values, follow the same idea.
# Text files with the final SFSs for the whole exome and the other four cases (see below), are available in the "data" directory of this repository.

for i in AW BC LB PG
do
	echo $i
	case $i in
		AW)
		projection=26
		;;
		BC)
		projection=12
		;;
		LB)
		projection=18
		;;
		PG)
		projection=28
		;;
	esac
	dadi-cli GenerateFs --vcf $VCFDIR/$i.S.vcf.gz  --pop-info $VCFDIR/${i}pop --pop-ids $i --projections $projection --output $SFSDIR/$i.S.exome.sfs
	dadi-cli GenerateFs --vcf $VCFDIR/$i.NS.vcf.gz --pop-info $VCFDIR/${i}pop --pop-ids $i --projections $projection --output $SFSDIR/$i.NS.exome.sfs
done


###########################################################################################
######################## #1 --> HC, exome, Nmax, gamma&neugamma ###########################
###########################################################################################

cd $HERE

###############
# DEMOGRAPHY  #
###############

# 2 EP #

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.exome.sfs 'two_epoch' './1_HC_Nmax_exome/AW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.exome.sfs 'two_epoch' './1_HC_Nmax_exome/BC/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.exome.sfs 'two_epoch' './1_HC_Nmax_exome/LB/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.exome.sfs 'two_epoch' './1_HC_Nmax_exome/PG/demo_2ep'

# 3 EP #

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.exome.sfs 'three_epoch' './1_HC_Nmax_exome/AW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.exome.sfs 'three_epoch' './1_HC_Nmax_exome/BC/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.exome.sfs 'three_epoch' './1_HC_Nmax_exome/LB/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.exome.sfs 'three_epoch' './1_HC_Nmax_exome/PG/demo_3ep'

#############
#### DFE ####
#############

# First Step of DFE inference - Generate a demographics informed precomputed spectra for each species/population.

python3 DFE1D_refspectra.py two_epoch '0.1884,0.0137' 26 './1_HC_Nmax_exome/AW/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.3320,0.1067' 12 './1_HC_Nmax_exome/BC/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.1994,0.1860' 18 './1_HC_Nmax_exome/LB/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.0356,0.0815' 28 './1_HC_Nmax_exome/PG/cache2ep'
python3 DFE1D_refspectra.py three_epoch '0.0008,0.1079,0.0008,0.01968' 28 './1_HC_Nmax_exome/PG/cache3ep'

# Last step of DFE inference - Infer DFE of NS mutations

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './1_HC_Nmax_exome/AW/cache2ep_DFESpectrum.bpkl' 'gamma' '10380.35531' './1_HC_Nmax_exome/AW/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './1_HC_Nmax_exome/BC/cache2ep_DFESpectrum.bpkl' 'gamma' '7357.450526' './1_HC_Nmax_exome/BC/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './1_HC_Nmax_exome/LB/cache2ep_DFESpectrum.bpkl' 'gamma' '11847.98314' './1_HC_Nmax_exome/LB/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './1_HC_Nmax_exome/PG/cache2ep_DFESpectrum.bpkl' 'gamma' '32021.06728' './1_HC_Nmax_exome/PG/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './1_HC_Nmax_exome/PG/cache3ep_DFESpectrum.bpkl' 'gamma' '13505.36804' './1_HC_Nmax_exome/PG/dfe/gamma_3ep'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './1_HC_Nmax_exome/AW/cache2ep_DFESpectrum.bpkl' 'neugamma' '10380.35531' './1_HC_Nmax_exome/AW/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './1_HC_Nmax_exome/BC/cache2ep_DFESpectrum.bpkl' 'neugamma' '7357.450526' './1_HC_Nmax_exome/BC/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './1_HC_Nmax_exome/LB/cache2ep_DFESpectrum.bpkl' 'neugamma' '11847.98314' './1_HC_Nmax_exome/LB/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 200 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './1_HC_Nmax_exome/PG/cache2ep_DFESpectrum.bpkl' 'neugamma' '32021.06728' './1_HC_Nmax_exome/PG/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './1_HC_Nmax_exome/PG/cache3ep_DFESpectrum.bpkl' 'neugamma' '13505.36804' './1_HC_Nmax_exome/PG/dfe/neugamma_3ep'


##########################################################################################
############################# #2 --> HC, exome, Neq, g+ng ################################
##########################################################################################

cd $HERE

##############
# DEMOGRAPHY #
##############

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.Neq.sfs 'two_epoch' './2_HC_Neq/AW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.Neq.sfs 'two_epoch' './2_HC_Neq/BC/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.Neq.sfs 'two_epoch' './2_HC_Neq/LB/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.Neq.sfs 'two_epoch' './2_HC_Neq/PG/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.Neq.sfs 'three_epoch' './2_HC_Neq/AW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.Neq.sfs 'three_epoch' './2_HC_Neq/BC/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.Neq.sfs 'three_epoch' './2_HC_Neq/LB/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.Neq.sfs 'three_epoch' './2_HC_Neq/PG/demo_3ep'

#############
#### DFE ####
#############

# First Step of DFE inference - Generate a demographics informed precomputed spectra for each species/population.

cd $HERE

python3 DFE1D_refspectra.py two_epoch '0.0239,0.0013' 12 './2_HC_Neq/AW/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.3320,0.1067' 12 './2_HC_Neq/BC/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.2295,0.1579' 12 './2_HC_Neq/LB/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.0411,0.0822' 12 './2_HC_Neq/PG/cache2ep'
python3 DFE1D_refspectra.py three_epoch '0.0033,0.1284,0.0036,0.0161' 12 './2_HC_Neq/PG/cache3ep'

# Last step of DFE inference - Infer DFE of NS mutations

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.Neq.sfs './2_HC_Neq/AW/cache2ep_DFESpectrum.bpkl' 'gamma' '10456.71859' './2_HC_Neq/AW/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.Neq.sfs './2_HC_Neq/BC/cache2ep_DFESpectrum.bpkl' 'gamma' '7357.744828' './2_HC_Neq/BC/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.Neq.sfs './2_HC_Neq/LB/cache2ep_DFESpectrum.bpkl' 'gamma' '10236.46092' './2_HC_Neq/LB/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.Neq.sfs './2_HC_Neq/PG/cache2ep_DFESpectrum.bpkl' 'gamma' '25469.65981' './2_HC_Neq/PG/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.Neq.sfs './2_HC_Neq/PG/cache3ep_DFESpectrum.bpkl' 'gamma' '13632.73012' './2_HC_Neq/PG/dfe/gamma_3ep'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.Neq.sfs './2_HC_Neq/AW/cache2ep_DFESpectrum.bpkl' 'neugamma' '10456.71859' './2_HC_Neq/AW/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.Neq.sfs './2_HC_Neq/BC/cache2ep_DFESpectrum.bpkl' 'neugamma' '7357.744828' './2_HC_Neq/BC/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.Neq.sfs './2_HC_Neq/LB/cache2ep_DFESpectrum.bpkl' 'neugamma' '10236.46092' './2_HC_Neq/LB/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.Neq.sfs './2_HC_Neq/PG/cache2ep_DFESpectrum.bpkl' 'neugamma' '25469.65981' './2_HC_Neq/PG/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.Neq.sfs './2_HC_Neq/PG/cache3ep_DFESpectrum.bpkl' 'neugamma' '13632.73012' './2_HC_Neq/PG/dfe/neugamma_3ep'

##########################################################################################
############################# #3 --> HC, gene sets, g+ng #################################
##########################################################################################

cd $HERE

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '3040579' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.immune.sfs 'two_epoch' './3_genesets/AW/immune/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '2978994' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.immune.sfs 'two_epoch' './3_genesets/BC/immune/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '2895641' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.immune.sfs 'two_epoch' './3_genesets/LB/immune/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '3015569' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.immune.sfs 'two_epoch' './3_genesets/PG/immune/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '3040579' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.immune.sfs 'three_epoch' './3_genesets/AW/immune/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '2978994' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.immune.sfs 'three_epoch' './3_genesets/BC/immune/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '2895641' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.immune.sfs 'three_epoch' './3_genesets/LB/immune/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '3015569' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.immune.sfs 'three_epoch' './3_genesets/PG/immune/demo_3ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '4729731' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.neuro.sfs 'two_epoch' './3_genesets/AW/neuro/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '4633933' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.neuro.sfs 'two_epoch' './3_genesets/BC/neuro/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '4504275' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.neuro.sfs 'two_epoch' './3_genesets/LB/neuro/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '4690828' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.neuro.sfs 'two_epoch' './3_genesets/PG/neuro/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '4729731' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.neuro.sfs 'three_epoch' './3_genesets/AW/neuro/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '4633933' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.neuro.sfs 'three_epoch' './3_genesets/BC/neuro/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '4504275' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.neuro.sfs 'three_epoch' './3_genesets/LB/neuro/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '4690828' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.neuro.sfs 'three_epoch' './3_genesets/PG/neuro/demo_3ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '8606538' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.domest.sfs 'two_epoch' './3_genesets/AW/domest/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '8432218' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.domest.sfs 'two_epoch' './3_genesets/BC/domest/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '8196284' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.domest.sfs 'two_epoch' './3_genesets/LB/domest/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '8535747' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.domest.sfs 'two_epoch' './3_genesets/PG/domest/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '8606538' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.domest.sfs 'three_epoch' './3_genesets/AW/domest/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '8432218' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.domest.sfs 'three_epoch' './3_genesets/BC/domest/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '8196284' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.domest.sfs 'three_epoch' './3_genesets/LB/domest/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '8535747' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.domest.sfs 'three_epoch' './3_genesets/PG/domest/demo_3ep'


#############
#### DFE ####
#############

# First Step of DFE inference - Generate a demographics informed precomputed spectra for each species/population.

cd $HERE

python3 DFE1D_refspectra.py two_epoch '0.29851376,0.030063324' 26 './3_genesets/AW/immune/cache2ep'
python3 DFE1D_refspectra.py two_epoch '9.86E-05,3.52E-06' 26 './3_genesets/AW/neuro/cache2ep'
python3 DFE1D_refspectra.py three_epoch '4.120700773,0.734370686,1.727464941,0.184651193' 26 './3_genesets/AW/domest/cache3ep'
python3 DFE1D_refspectra.py two_epoch '0.170690473,0.029488129' 12 './3_genesets/BC/immune/cache2ep'
python3 DFE1D_refspectra.py three_epoch '0.008654746,3.387762787,0.011631614,0.013657794' 12 './3_genesets/BC/neuro/cache3ep'
python3 DFE1D_refspectra.py two_epoch '0.112672094,0.209952771' 12 './3_genesets/BC/domest/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.300321909,0.17312621' 18 './3_genesets/LB/immune/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.311823054,0.198585014' 18 './3_genesets/LB/neuro/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.314037126,0.123708776' 18 './3_genesets/LB/domest/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.070612287,0.120820319' 24 './3_genesets/PG/immune/cache2ep'
python3 DFE1D_refspectra.py three_epoch '0.009308742,0.279635331,0.008443861,0.012645507' 24 './3_genesets/PG/neuro/cache3ep'
python3 DFE1D_refspectra.py three_epoch '0.001836578,0.162209946,0.001417976,0.020704445' 24 './3_genesets/PG/domest/cache3ep'

# Last step of DFE inference - Infer DFE of NS mutations

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '3040579' --NS_S_scaling '2.21' $SFSDIR/AW.NS.immune.sfs './3_genesets/AW/immune/cache2ep_DFESpectrum.bpkl' 'gamma' '610.5040' './3_genesets/AW/immune/dfe/gamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '4729731' --NS_S_scaling '2.21' $SFSDIR/AW.NS.neuro.sfs  './3_genesets/AW/neuro/cache2ep_DFESpectrum.bpkl' 'gamma' '793.4689' './3_genesets/AW/neuro/dfe/gamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '8606538' --NS_S_scaling '2.21' $SFSDIR/AW.NS.domest.sfs './3_genesets/AW/domest/cache3ep_DFESpectrum.bpkl' 'gamma' '838.4166' './3_genesets/AW/domest/dfe/gamma_domest'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '3040579' --NS_S_scaling '2.21' $SFSDIR/AW.NS.immune.sfs './3_genesets/AW/immune/cache2ep_DFESpectrum.bpkl' 'neugamma' '610.5040' './3_genesets/AW/immune/dfe/neugamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '4729731' --NS_S_scaling '2.21' $SFSDIR/AW.NS.neuro.sfs  './3_genesets/AW/neuro/cache2ep_DFESpectrum.bpkl' 'neugamma' '793.4689' './3_genesets/AW/neuro/dfe/neugamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '8606538' --NS_S_scaling '2.21' $SFSDIR/AW.NS.domest.sfs './3_genesets/AW/domest/cache3ep_DFESpectrum.bpkl' 'neugamma' '838.4166' './3_genesets/AW/domest/dfe/neugamma_domest'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '2978994' --NS_S_scaling '2.21' $SFSDIR/BC.NS.immune.sfs './3_genesets/BC/immune/cache2ep_DFESpectrum.bpkl' 'gamma' '399.3353' './3_genesets/BC/immune/dfe/gamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '4633933' --NS_S_scaling '2.21' $SFSDIR/BC.NS.neuro.sfs  './3_genesets/BC/neuro/cache3ep_DFESpectrum.bpkl' 'gamma' '1695.1118' './3_genesets/BC/neuro/dfe/gamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '8432218' --NS_S_scaling '2.21' $SFSDIR/BC.NS.domest.sfs './3_genesets/BC/domest/cache2ep_DFESpectrum.bpkl' 'gamma' '3702.2009' './3_genesets/BC/domest/dfe/gamma_domest'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '2978994' --NS_S_scaling '2.21' $SFSDIR/BC.NS.immune.sfs './3_genesets/BC/immune/cache2ep_DFESpectrum.bpkl' 'neugamma' '399.3353' './3_genesets/BC/immune/dfe/neugamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '4633933' --NS_S_scaling '2.21' $SFSDIR/BC.NS.neuro.sfs  './3_genesets/BC/neuro/cache3ep_DFESpectrum.bpkl' 'neugamma' '1695.1118' './3_genesets/BC/neuro/dfe/neugamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '8432218' --NS_S_scaling '2.21' $SFSDIR/BC.NS.domest.sfs './3_genesets/BC/domest/cache2ep_DFESpectrum.bpkl' 'neugamma' '3702.2009' './3_genesets/BC/domest/dfe/neugamma_domest'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '2895641' --NS_S_scaling '2.21' $SFSDIR/LB.NS.immune.sfs './3_genesets/LB/immune/cache2ep_DFESpectrum.bpkl' 'gamma' '525.2897' './3_genesets/LB/immune/dfe/gamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '4504275' --NS_S_scaling '2.21' $SFSDIR/LB.NS.neuro.sfs  './3_genesets/LB/neuro/cache2ep_DFESpectrum.bpkl' 'gamma' '715.5249' './3_genesets/LB/neuro/dfe/gamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '8196284' --NS_S_scaling '2.21' $SFSDIR/LB.NS.domest.sfs './3_genesets/LB/domest/cache2ep_DFESpectrum.bpkl' 'gamma' '1223.1277' './3_genesets/LB/domest/dfe/gamma_domest'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '2895641' --NS_S_scaling '2.21' $SFSDIR/LB.NS.immune.sfs './3_genesets/LB/immune/cache2ep_DFESpectrum.bpkl' 'neugamma' '525.2897' './3_genesets/LB/immune/dfe/neugamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '4504275' --NS_S_scaling '2.21' $SFSDIR/LB.NS.neuro.sfs  './3_genesets/LB/neuro/cache2ep_DFESpectrum.bpkl' 'neugamma' '715.5249' './3_genesets/LB/neuro/dfe/neugamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '8196284' --NS_S_scaling '2.21' $SFSDIR/LB.NS.domest.sfs './3_genesets/LB/domest/cache2ep_DFESpectrum.bpkl' 'neugamma' '1223.1277' './3_genesets/LB/domest/dfe/neugamma_domest'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '3015569' --NS_S_scaling '2.21' $SFSDIR/PG.NS.immune.sfs './3_genesets/PG/immune/cache2ep_DFESpectrum.bpkl' 'gamma' '1007.8588' './3_genesets/PG/immune/dfe/gamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '4690828' --NS_S_scaling '2.21' $SFSDIR/PG.NS.neuro.sfs  './3_genesets/PG/neuro/cache3ep_DFESpectrum.bpkl' 'gamma' '814.5026' './3_genesets/PG/neuro/dfe/gamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '8535747' --NS_S_scaling '2.21' $SFSDIR/PG.NS.domest.sfs './3_genesets/PG/domest/cache3ep_DFESpectrum.bpkl' 'gamma' '1541.7875' './3_genesets/PG/domest/dfe/gamma_domest'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '3015569' --NS_S_scaling '2.21' $SFSDIR/PG.NS.immune.sfs './3_genesets/PG/immune/cache2ep_DFESpectrum.bpkl' 'neugamma' '1007.8588' './3_genesets/PG/immune/dfe/neugamma_immune'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '4690828' --NS_S_scaling '2.21' $SFSDIR/PG.NS.neuro.sfs  './3_genesets/PG/neuro/cache3ep_DFESpectrum.bpkl' 'neugamma' '814.5026' './3_genesets/PG/neuro/dfe/neugamma_neuro'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '8535747' --NS_S_scaling '2.21' $SFSDIR/PG.NS.domest.sfs './3_genesets/PG/domest/cache3ep_DFESpectrum.bpkl' 'neugamma' '1541.7875' './3_genesets/PG/domest/dfe/neugamma_domest'


##########################################################################################
############################# #4 --> Low coverage samples ################################
##########################################################################################

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'MD' --mu '5.625E-09' --Lcds '15136299' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/MD.S.exome.sfs 'two_epoch' './4_LC/MD/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'MW' --mu '5.625E-09' --Lcds '15860106' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/MW.S.exome.sfs 'two_epoch' './4_LC/MW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'TM' --mu '5.625E-09' --Lcds '10852734' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/TM.S.exome.sfs 'two_epoch' './4_LC/TM/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'MD' --mu '5.625E-09' --Lcds '15136299' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/MD.S.exome.sfs 'three_epoch' './4_LC/MD/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'MW' --mu '5.625E-09' --Lcds '15860106' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/MW.S.exome.sfs 'three_epoch' './4_LC/MW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'TM' --mu '5.625E-09' --Lcds '10852734' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/TM.S.exome.sfs 'three_epoch' './4_LC/TM/demo_3ep'

#############
#### DFE ####
#############

# First Step of DFE inference - Generate a demographics informed precomputed spectra for each species/population.

python3 DFE1D_refspectra.py three_epoch '0.495852988,0.001420376,0.918612047,1.90E-05' 32 './4_LC/MD/cache3ep'
python3 DFE1D_refspectra.py two_epoch '3.995219504,3.766764239' 16 './4_LC/MW/cache2ep'
python3 DFE1D_refspectra.py two_epoch '0.432626904,0.073902898' 14 './4_LC/TM/cache2ep'

# Last step of DFE inference - Infer DFE of NS mutations

OUT=/Users/admin/Desktop/K9DFE_Unil_local/RUNS/4_LC
SFSDIR=/Users/admin/Desktop/K9DFE_Unil_local/DATA

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'MD' --mu '5.625E-09' --Lcds '15136299' --NS_S_scaling '2.21' $SFSDIR/MD.NS.exome.sfs './4_LC/MD/cache3ep_DFESpectrum.bpkl' 'gamma' '9037.89268' './4_LC/MD/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'MW' --mu '5.625E-09' --Lcds '15860106' --NS_S_scaling '2.21' $SFSDIR/MW.NS.exome.sfs './4_LC/MW/cache2ep_DFESpectrum.bpkl' 'gamma' '2594.52022' './4_LC/MW/dfe/gamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'TM' --mu '5.625E-09' --Lcds '10852734' --NS_S_scaling '2.21' $SFSDIR/TM.NS.exome.sfs './4_LC/TM/cache2ep_DFESpectrum.bpkl' 'gamma' '3917.51262' './4_LC/TM/dfe/gamma'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'MD' --mu '5.625E-09' --Lcds '15136299' --NS_S_scaling '2.21' $SFSDIR/MD.NS.exome.sfs './4_LC/MD/cache3ep_DFESpectrum.bpkl' 'neugamma' '9037.89268' './4_LC/MD/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'MW' --mu '5.625E-09' --Lcds '15860106' --NS_S_scaling '2.21' $SFSDIR/MW.NS.exome.sfs './4_LC/MW/cache2ep_DFESpectrum.bpkl' 'neugamma' '2594.52022' './4_LC/MW/dfe/neugamma'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'TM' --mu '5.625E-09' --Lcds '10852734' --NS_S_scaling '2.21' $SFSDIR/TM.NS.exome.sfs './4_LC/TM/cache2ep_DFESpectrum.bpkl' 'neugamma' '3917.51262' './4_LC/TM/dfe/neugamma'

##########################################################################################
########################## Runs #5 --> diff mutation and NSS #############################
##########################################################################################

# high mutation rate
# low mutation rate
# human NS:S ratio
# assuming 3-epoch for all four HC population samples (below)

###############
# DEMOGRAPHY  #
###############

# 2 EP #

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.31' --initval '0.3,0.1' $SFSDIR/AW.S.exome.sfs 'two_epoch' './5_other/humanNSS/AW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.31' --initval '0.3,0.1' $SFSDIR/BC.S.exome.sfs 'two_epoch' './5_other/humanNSS/BC/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.31' --initval '0.3,0.1' $SFSDIR/LB.S.exome.sfs 'two_epoch' './5_other/humanNSS/LB/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' --initval '0.3,0.1' $SFSDIR/PG.S.exome.sfs 'two_epoch' './5_other/humanNSS/PG/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '6.730e-9' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.exome.sfs 'two_epoch' './5_other/highmu/AW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '6.730e-9' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.exome.sfs 'two_epoch' './5_other/highmu/BC/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '6.730e-9' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.exome.sfs 'two_epoch' './5_other/highmu/LB/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.exome.sfs 'two_epoch' './5_other/highmu/PG/demo_2ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '3.000E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/AW.S.exome.sfs 'two_epoch' './5_other/lowmu/AW/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '3.000E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/BC.S.exome.sfs 'two_epoch' './5_other/lowmu/BC/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '3.000E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/LB.S.exome.sfs 'two_epoch' './5_other/lowmu/LB/demo_2ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '3.000E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.1' $SFSDIR/PG.S.exome.sfs 'two_epoch' './5_other/lowmu/PG/demo_2ep'

# 3 EP #

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.31' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.exome.sfs 'three_epoch' './5_other/humanNSS/AW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.31' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.exome.sfs 'three_epoch' './5_other/humanNSS/BC/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.31' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.exome.sfs 'three_epoch' './5_other/humanNSS/LB/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.exome.sfs 'three_epoch' './5_other/humanNSS/PG/demo_3ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '6.730e-9' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.exome.sfs 'three_epoch' './5_other/highmu/AW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '6.730e-9' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.exome.sfs 'three_epoch' './5_other/highmu/BC/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '6.730e-9' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.exome.sfs 'three_epoch' './5_other/highmu/LB/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.exome.sfs 'three_epoch' './5_other/highmu/PG/demo_3ep'

python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'AW' --mu '3.000E-09' --Lcds '20644027' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/AW.S.exome.sfs 'three_epoch' './5_other/lowmu/AW/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'BC' --mu '3.000E-09' --Lcds '20225895' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/BC.S.exome.sfs 'three_epoch' './5_other/lowmu/BC/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'LB' --mu '3.000E-09' --Lcds '19659972' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/LB.S.exome.sfs 'three_epoch' './5_other/lowmu/LB/demo_3ep'
python3 Demog1D_sizechangeFIM.py --Nrun 100 --pop 'PG' --mu '3.000E-09' --Lcds '20474224' --NS_S_scaling '2.21' --initval '0.3,0.3,0.1,0.1' $SFSDIR/PG.S.exome.sfs 'three_epoch' './5_other/lowmu/PG/demo_3ep'

#############
#### DFE ####
#############

# First Step - Generate a demographics informed precomputed spectra for each species/population.

### 2-epoch 

python3 DFE1D_refspectra.py two_epoch '0.187428927,0.013601214' 26 './5_other/highmu/AW/cache2ep_highmu'
python3 DFE1D_refspectra.py two_epoch '0.331896584,0.106628021' 12 './5_other/highmu/BC/cache2ep_highmu'
python3 DFE1D_refspectra.py two_epoch '0.199002809,0.186173734' 18 './5_other/highmu/LB/cache2ep_highmu'
python3 DFE1D_refspectra.py two_epoch '0.045192570,0.092101313' 28 './5_other/highmu/PG/cache2ep_highmu'
python3 DFE1D_refspectra.py two_epoch '0.188324237,0.013685284' 26 './5_other/humanNSS/AW/cache2ep_humanNSS'
python3 DFE1D_refspectra.py two_epoch '0.331926452,0.106655237' 12 './5_other/humanNSS/BC/cache2ep_humanNSS'
python3 DFE1D_refspectra.py two_epoch '0.199348823,0.186020281' 18 './5_other/humanNSS/LB/cache2ep_humanNSS'
python3 DFE1D_refspectra.py two_epoch '0.045587801,0.092480275' 28 './5_other/humanNSS/PG/cache2ep_humanNSS'
python3 DFE1D_refspectra.py two_epoch '0.187155178,0.013580639' 26 './5_other/lowmu/AW/cache2ep_lowmu'
python3 DFE1D_refspectra.py two_epoch '0.331794941,0.106519964' 12 './5_other/lowmu/BC/cache2ep_lowmu'
python3 DFE1D_refspectra.py two_epoch '0.199186973,0.186120156' 18 './5_other/lowmu/LB/cache2ep_lowmu'
python3 DFE1D_refspectra.py two_epoch '0.046142993,0.093044406' 28 './5_other/lowmu/PG/cache2ep_lowmu'
python3 DFE1D_refspectra.py three_epoch '0.000284203,0.075451311,0.000386971,0.01433693' 28 './5_other/highmu/PG/cache3ep_highmu'

### 3-epoch for PG

python3 DFE1D_refspectra.py three_epoch '0.000133003,0.098539087,0.00014717,0.018319254' 28 './5_other/humanNSS/PG/cache3ep_humanNSS'
python3 DFE1D_refspectra.py three_epoch '2.10E-05,0.09896017,2.27E-05,0.019180262' 28 './5_other/lowmu/PG/cache3ep_lowmu'

# Last step of DFE inference - Infer DFE of NS mutations

# Human NS:S

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.31' $SFSDIR/AW.NS.exome.sfs './5_other/humanNSS/AW/cache2ep_humanNSS_DFESpectrum.bpkl' 'gamma' '10379.68389' './5_other/humanNSS/AW/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.31' $SFSDIR/BC.NS.exome.sfs './5_other/humanNSS/BC/cache2ep_humanNSS_DFESpectrum.bpkl' 'gamma' '7357.148503' './5_other/humanNSS/BC/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.31' $SFSDIR/LB.NS.exome.sfs './5_other/humanNSS/LB/cache2ep_humanNSS_DFESpectrum.bpkl' 'gamma' '11850.11164' './5_other/humanNSS/LB/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' $SFSDIR/PG.NS.exome.sfs './5_other/humanNSS/PG/cache2ep_humanNSS_DFESpectrum.bpkl' 'gamma' '24950.14526' './5_other/humanNSS/PG/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' $SFSDIR/PG.NS.exome.sfs './5_other/humanNSS/PG/cache3ep_humanNSS_DFESpectrum.bpkl' 'gamma' '14711.51874' './5_other/humanNSS/PG/dfe/gamma_3ep'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.31' $SFSDIR/AW.NS.exome.sfs './5_other/humanNSS/AW/cache2ep_humanNSS_DFESpectrum.bpkl' 'neugamma' '10379.68389' './5_other/humanNSS/AW/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.31' $SFSDIR/BC.NS.exome.sfs './5_other/humanNSS/BC/cache2ep_humanNSS_DFESpectrum.bpkl' 'neugamma' '7357.148503' './5_other/humanNSS/BC/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.31' $SFSDIR/LB.NS.exome.sfs './5_other/humanNSS/LB/cache2ep_humanNSS_DFESpectrum.bpkl' 'neugamma' '11850.11164' './5_other/humanNSS/LB/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' $SFSDIR/PG.NS.exome.sfs './5_other/humanNSS/PG/cache2ep_humanNSS_DFESpectrum.bpkl' 'neugamma' '24950.14526' './5_other/humanNSS/PG/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.31' $SFSDIR/PG.NS.exome.sfs './5_other/humanNSS/PG/cache3ep_humanNSS_DFESpectrum.bpkl' 'neugamma' '14711.51874' './5_other/humanNSS/PG/dfe/neugamma_3ep'

# high mu

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '6.730e-9' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './5_other/highmu/AW/cache2ep_highmu_DFESpectrum.bpkl' 'gamma' '10379.7754' './5_other/highmu/AW/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '6.730e-9' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './5_other/highmu/BC/cache2ep_highmu_DFESpectrum.bpkl' 'gamma' '7357.00131' './5_other/highmu/BC/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '6.730e-9' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './5_other/highmu/LB/cache2ep_highmu_DFESpectrum.bpkl' 'gamma' '11872.3231' './5_other/highmu/LB/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/highmu/PG/cache2ep_highmu_DFESpectrum.bpkl' 'gamma' '25171.9793' './5_other/highmu/PG/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/highmu/PG/cache3ep_highmu_DFESpectrum.bpkl' 'gamma' '19037.7512' './5_other/highmu/PG/dfe/gamma_3ep'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '6.730e-9' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './5_other/highmu/AW/cache2ep_highmu_DFESpectrum.bpkl' 'neugamma' '10379.7754' './5_other/highmu/AW/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '6.730e-9' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './5_other/highmu/BC/cache2ep_highmu_DFESpectrum.bpkl' 'neugamma' '7357.00131' './5_other/highmu/BC/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '6.730e-9' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './5_other/highmu/LB/cache2ep_highmu_DFESpectrum.bpkl' 'neugamma' '11872.3231' './5_other/highmu/LB/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/highmu/PG/cache2ep_highmu_DFESpectrum.bpkl' 'neugamma' '25171.9793' './5_other/highmu/PG/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '6.730e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/highmu/PG/cache3ep_highmu_DFESpectrum.bpkl' 'neugamma' '19037.7512' './5_other/highmu/PG/dfe/neugamma_3ep'

# low mu

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '3.000e-9' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './5_other/lowmu/AW/cache2ep_lowmu_DFESpectrum.bpkl' 'gamma' '10380.2856' './5_other/lowmu/AW/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '3.000e-9' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './5_other/lowmu/BC/cache2ep_lowmu_DFESpectrum.bpkl' 'gamma' '7356.27602' './5_other/lowmu/BC/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '3.000e-9' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './5_other/lowmu/LB/cache2ep_lowmu_DFESpectrum.bpkl' 'gamma' '11861.4255' './5_other/lowmu/LB/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '3.000e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/lowmu/PG/cache2ep_lowmu_DFESpectrum.bpkl' 'gamma' '24658.5916' './5_other/lowmu/PG/dfe/gamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '3.000e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/lowmu/PG/cache3ep_lowmu_DFESpectrum.bpkl' 'gamma' '14711.5187' './5_other/lowmu/PG/dfe/gamma_3ep'

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '3.000e-9' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './5_other/lowmu/AW/cache2ep_lowmu_DFESpectrum.bpkl' 'neugamma' '10380.2856' './5_other/lowmu/AW/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '3.000e-9' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './5_other/lowmu/BC/cache2ep_lowmu_DFESpectrum.bpkl' 'neugamma' '7356.27602' './5_other/lowmu/BC/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '3.000e-9' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './5_other/lowmu/LB/cache2ep_lowmu_DFESpectrum.bpkl' 'neugamma' '11861.4255' './5_other/lowmu/LB/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '3.000e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/lowmu/PG/cache2ep_lowmu_DFESpectrum.bpkl' 'neugamma' '24658.5916' './5_other/lowmu/PG/dfe/neugamma_2ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '3.000e-9' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './5_other/lowmu/PG/cache3ep_lowmu_DFESpectrum.bpkl' 'neugamma' '14711.5187' './5_other/lowmu/PG/dfe/neugamma_3ep'

#### Gamma DFE considering 3ep demography for all 4 HC pops considering same dataset as #1.

#############
#### DFE ####
#############

# First Step - Generate a demographics informed precomputed spectra for each species/population.

python3 DFE1D_refspectra.py three_epoch '2.0805,0.4944,3.0084,0.0570' 26 './1_HC_Nmax_exome/AW/cache3ep'
python3 DFE1D_refspectra.py three_epoch '0.2342,0.00004,0.3809,0.000001' 12 './1_HC_Nmax_exome/BC/cache3ep'
python3 DFE1D_refspectra.py three_epoch '0.1588,0.1506,0.1920,0.0062' 18 './1_HC_Nmax_exome/LB/cache3ep'
python3 DFE1D_refspectra.py three_epoch '0.0008,0.1079,0.0008,0.0196' 28 './1_HC_Nmax_exome/PG/cache3ep'

# Last step of DFE inference - Infer DFE of NS mutations

python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'AW' --mu '5.625E-09' --Lcds '20644027' --NS_S_scaling '2.21' $SFSDIR/AW.NS.exome.sfs './1_HC_Nmax_exome/AW/cache3ep_DFESpectrum.bpkl' 'gamma' '5818.269478' './1_HC_Nmax_exome/AW/dfe/gamma_3ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'BC' --mu '5.625E-09' --Lcds '20225895' --NS_S_scaling '2.21' $SFSDIR/BC.NS.exome.sfs './1_HC_Nmax_exome/BC/cache3ep_DFESpectrum.bpkl' 'gamma' '16030.39622' './1_HC_Nmax_exome/BC/dfe/gamma_3ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'LB' --mu '5.625E-09' --Lcds '19659972' --NS_S_scaling '2.21' $SFSDIR/LB.NS.exome.sfs './1_HC_Nmax_exome/LB/cache3ep_DFESpectrum.bpkl' 'gamma' '15250.47593' './1_HC_Nmax_exome/LB/dfe/gamma_3ep'
python3 DFE1D_inferenceFIM.py --Nrun 50 --pop 'PG' --mu '5.625E-09' --Lcds '20474224' --NS_S_scaling '2.21' $SFSDIR/PG.NS.exome.sfs './1_HC_Nmax_exome/PG/cache3ep_DFESpectrum.bpkl' 'gamma' '13505.36804' './1_HC_Nmax_exome/PG/dfe/gamma_3ep'
