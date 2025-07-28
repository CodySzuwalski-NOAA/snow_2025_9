## GMACS Version 2.20.22; ** AEP & WTS **; Compiled 2025-03-31

# Block structure
# Number of block groups
2
# Block structure (number of blocks per block group)
3 # block group 1
42 # block group 2
# The blocks
#Block 1: 
2018 2018 # block_group_1_block_1
2019 2019 # block_group_1_block_2
2020 2020 # block_group_1_block_3
#Block 2: 
1983 1983 # block_group_2_block_1
1984 1984 # block_group_2_block_2
1985 1985 # block_group_2_block_3
1986 1986 # block_group_2_block_4
1987 1987 # block_group_2_block_5
1988 1988 # block_group_2_block_6
1989 1989 # block_group_2_block_7
1990 1990 # block_group_2_block_8
1991 1991 # block_group_2_block_9
1992 1992 # block_group_2_block_10
1993 1993 # block_group_2_block_11
1994 1994 # block_group_2_block_12
1995 1995 # block_group_2_block_13
1996 1996 # block_group_2_block_14
1997 1997 # block_group_2_block_15
1998 1998 # block_group_2_block_16
1999 1999 # block_group_2_block_17
2000 2000 # block_group_2_block_18
2001 2001 # block_group_2_block_19
2002 2002 # block_group_2_block_20
2003 2003 # block_group_2_block_21
2004 2004 # block_group_2_block_22
2005 2005 # block_group_2_block_23
2006 2006 # block_group_2_block_24
2007 2007 # block_group_2_block_25
2008 2008 # block_group_2_block_26
2009 2009 # block_group_2_block_27
2010 2010 # block_group_2_block_28
2011 2011 # block_group_2_block_29
2012 2012 # block_group_2_block_30
2013 2013 # block_group_2_block_31
2014 2014 # block_group_2_block_32
2015 2015 # block_group_2_block_33
2016 2016 # block_group_2_block_34
2017 2017 # block_group_2_block_35
2018 2018 # block_group_2_block_36
2019 2019 # block_group_2_block_37
2020 2020 # block_group_2_block_38
2021 2021 # block_group_2_block_39
2022 2022 # block_group_2_block_40
2023 2023 # block_group_2_block_41
2024 2024 # block_group_2_block_42

##  ------------------------------------------------------------------------------------ ##
##  OTHER  CONTROLS
##  ------------------------------------------------------------------------------------ ##
1982 # First year of recruitment estimation
2023 # Last year of recruitment estimation
   1 # Consider terminal molting (0 = off, 1 = on). If on, the calc_stock_recruitment_relationship() isn't called in the procedure
   1 # Phase for recruitment estimation
   2 # Phase for recruitment sex-ratio estimation
0.50 # Initial value for recruitment sex-ratio
   2 # Initial conditions (0 = Unfished, 1 = Steady-state fished, 2 = Free parameters, 3 = Free parameters (revised))
   1 # Reference size-class for initial conditons = 3
1.00 # Lambda (proportion of mature male biomass for SPR reference points)
   0 # Stock-Recruit-Relationship (0 = none, 1 = Beverton-Holt)
   0 # Use years specified to computed average sex ratio in the calculation of average recruitment for reference points (0 = off -i.e. Rec based on End year, 1 = on)
 200 # Years to compute equilibria
   5 # Phase for deviation parameters
1940 # First year of bias-correction
1950 # First full bias-correction
2050 # Last full bias-correction
2051 # Last year of bias-correction

# Expecting 98 theta parameters
# Core parameters
## Initial: Initial value for the parameter (must lie between lower and upper)
## Lower & Upper: Range for the parameter
## Phase: Set equal to a negative number not to estimate
## Prior type:
## 0: Uniform   - parameters are the range of the uniform prior
## 1: Normal    - parameters are the mean and sd
## 2: Lognormal - parameters are the mean and sd of the log
## 3: Beta      - parameters are the two beta parameters [see dbeta]
## 4: Gamma     - parameters are the two gamma parameters [see dgamma]
# Initial_value    Lower_bound    Upper_bound Phase Prior_type        Prior_1        Prior_2
    16.50000000   -10.00000000    20.00000000    -2          0   -10.00000000    20.00000000 # Log(R0)
    15.00000000   -10.00000000    30.00000000    -1          0    10.00000000    20.00000000 # Log(Rinitial)
    13.26245375   -10.00000000    30.00000000     1          0    10.00000000    20.00000000 # Log(Rbar)
    32.50000000     7.50000000    42.50000000    -4          0    32.50000000     2.25000000 # Recruitment_ra-males
     1.00000000     0.10000000    10.00000000    -4          0     0.10000000     5.00000000 # Recruitment_rb-males
     0.00000000   -10.00000000    10.00000000    -4          0     0.00000000    20.00000000 # Recruitment_ra-females
     0.00000000   -10.00000000    10.00000000    -3          0     0.00000000    20.00000000 # Recruitment_rb-females
    -0.90000000   -10.00000000     0.75000000    -4          0   -10.00000000     0.75000000 # log(SigmaR)
     0.75000000     0.20000000     1.00000000    -2          3     3.00000000     2.00000000 # Steepness
     0.01000000     0.00010000     1.00000000    -3          3     1.01000000     1.01000000 # Rho
     0.00000000   -10.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_1
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_2
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_3
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_4
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_5
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_6
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_7
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_8
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_9
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_10
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_11
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_12
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_13
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_14
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_15
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_16
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_17
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_18
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_19
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_20
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_21
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_mature_newshell_class_22
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_1
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_2
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_3
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_4
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_5
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_6
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_7
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_8
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_9
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_10
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_11
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_12
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_13
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_14
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_15
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_16
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_17
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_18
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_19
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_20
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_21
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_male_mature_immature_newshell_class_22
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_1
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_2
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_3
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_4
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_5
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_6
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_7
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_8
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_9
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_10
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_11
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_12
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_13
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_14
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_15
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_16
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_17
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_18
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_19
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_20
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_21
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_mature_newshell_class_22
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_1
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_2
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_3
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_4
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_5
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_6
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_7
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_8
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_9
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_10
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_11
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_12
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_13
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_14
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_15
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_16
     0.00000000   -20.00000000    25.00000000     1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_17
   -19.00000000   -20.00000000    25.00000000    -1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_18
   -19.00000000   -20.00000000    25.00000000    -1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_19
   -19.00000000   -20.00000000    25.00000000    -1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_20
   -19.00000000   -20.00000000    25.00000000    -1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_21
   -19.00000000   -20.00000000    25.00000000    -1          0    10.00000000    20.00000000 # Initial_logN_for_sex_female_mature_immature_newshell_class_22

##  ------------------------------------------------------------------------------------ ##
## Allometry
##  ------------------------------------------------------------------------------------ ##
# weight-at-length input  method:
## 1 = allometry  [w_l = a*l^b]
## 2 = vector by sex
## 3 = matrix by sex
2  #--selected method
     0.00000766     0.00001290     0.00002000     0.00002950     0.00004170     0.00005680     0.00007530     0.00009745     0.00012369     0.00015433     0.00018974     0.00023028     0.00027631     0.00032821     0.00038633     0.00045106     0.00052275     0.00060180     0.00068856     0.00078342     0.00088677     0.00099897
     0.00000766     0.00001290     0.00002000     0.00002950     0.00004170     0.00005680     0.00007530     0.00009745     0.00012369     0.00015433     0.00018974     0.00023028     0.00027631     0.00032821     0.00038633     0.00045106     0.00052275     0.00060180     0.00068856     0.00078342     0.00088677     0.00099897
     0.00000917     0.00001440     0.00002130     0.00002980     0.00004030     0.00005290     0.00006770     0.00008480     0.00010445     0.00012676     0.00015186     0.00017988     0.00021096     0.00024523     0.00028282     0.00032385     0.00036845     0.00041673     0.00046883     0.00052485     0.00058492     0.00064916
     0.00000917     0.00001440     0.00002130     0.00002980     0.00004030     0.00005290     0.00006770     0.00008480     0.00010445     0.00012676     0.00015186     0.00017988     0.00021096     0.00024523     0.00028282     0.00032385     0.00036845     0.00041673     0.00046883     0.00052485     0.00058492     0.00064916
##  ------------------------------------------------------------------------------------ ##

##  ------------------------------------------------------------------------------------ ##
## Proportion mature by sex and size
 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
##  ------------------------------------------------------------------------------------ ##

##  ------------------------------------------------------------------------------------ ##
# Proportion legal by sex and size
 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
##  ------------------------------------------------------------------------------------ ##


## ==================================================================================== ##
## GROWTH PARAMETER CONTROLS                                                            ##
## ==================================================================================== ##
## 
# Maximum number of size-classes to which recruitment must occur
 3 3
# Use functional maturity for terminally molting animals (0=no; 1=Yes)?
0
# Growth transition
##Type_1: Options for the growth matrix
#  1: Pre-specified growth transition matrix (requires molt probability)
#  2: Pre-specified size transition matrix (molt probability is ignored)
#  3: Growth increment is gamma distributed (requires molt probability)
#  4: Post-molt size is gamma distributed (requires molt probability)
#  5: Von Bert.: kappa varies among individuals (requires molt probability)
#  6: Von Bert.: Linf varies among individuals (requires molt probability)
#  7: Von Bert.: kappa and Linf varies among individuals (requires molt probability)
#  8: Growth increment is normally distributed (requires molt probability)
## Type_2: Options for the growth increment model matrix
#  1: Linear
##     1. intercept 
##     2. slope
##     3. gamma distribution scale parameter (on ln-scale)
#  2: Individual
#  3: Individual (Same as 2)
#  4: Power law for mean post-molt size (3 parameters)
##     1. ln-scale intercept 
##     2. ln-scale slope
##     3. gamma distribution scale parameter (on ln-scale)
#  5: Alternative power law for mean post-molt size (5 parameters)
##     1. reference (small) pre-molt size (constant: must have phase < 0)
##     2. mean post-molt size corresponding to 1.
##     3. reference (large) pre-molt size (constant: must have phase < 0)
##     4. mean post-molt size corresponding to 3.
##     5. gamma distribution scale parameter (on arithmetic scale)
#  Block: Block number for time-varying growth   
## Type_1 Type_2  Block
        4      1      0 
        4      1      0 
# Molt probability
# Type: Options for the molt probability function
#  0: Pre-specified
#  1: Constant at 1
#  2: Logistic
#  3: Individual
#  Block: Block number for time-varying growth   
## Type Block
      1     0 
      1     0 
# Mature probability
# Type: Options for the mature probability function
#  0: Pre-specified
#  1: Constant at 1
#  2: Logistic
#  3: Individual
# Block: Block number for time-varying growth   
## Type Block
      0     2 
      0     2 

## General parameter specificiations 
##  Initial: Initial value for the parameter (must lie between lower and upper)
##  Lower & Upper: Range for the parameter
##  Prior type:
##   0: Uniform   - parameters are the range of the uniform prior
##   1: Normal    - parameters are the mean and sd
##   2: Lognormal - parameters are the mean and sd of the log
##   3: Beta      - parameters are the two beta parameters [see dbeta]
##   4: Gamma     - parameters are the two gamma parameters [see dgamma]
##  Phase: Set equal to a negative number not to estimate
##  Relative: 0: absolute; 1 relative 
##  Block: Block number for time-varying selectivity   
##  Block_fn: 0:absolute values; 1:exponential
##  Env_L: Environmental link - options are 0:none; 1:additive; 2:multiplicative; 3:exponential
##  EnvL_var: Environmental variable
##  RW: 0 for no random walk changes; 1 otherwise
##  RW_blk: Block number for random walks
##  Sigma_RW: Sigma used for the random walk

### Parameter inputs for growth transitions
# Inputs for sex * type 1
# MAIN PARS: Initial  Lower_bound  Upper_bound Prior_type       Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
             2.049000    -5.000000    20.000000          0     2.049000     1.000000      3      0      0      0      0      0      0   0.3000 # Alpha_base_male
            -0.225800    -1.000000     0.000000          0    -0.225800     0.500000      3      0      0      0      0      0      0   0.3000 # Beta_base_male
             0.250000     0.001000     5.000000          0     0.000000   999.000000     -3      0      0      0      0      0      0   0.3000 # Gscale_base_male
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

### Parameter inputs for growth transitions
# Inputs for sex * type 2
# MAIN PARS: Initial  Lower_bound  Upper_bound Prior_type       Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
            -1.153900    -5.000000    10.000000          0    -1.153900     1.000000      3      0      0      0      0      0      0   0.3000 # Alpha_base_female
            -0.338900    -1.000000     0.000000          0    -0.338900     0.500000      3      0      0      0      0      0      0   0.3000 # Beta_base_female
             0.250000     0.001000     5.000000          0     0.000000   999.000000     -3      0      0      0      0      0      0   0.3000 # Gscale_base_female
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

### Parameter inputs for probability of molting
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

### Parameter inputs for probability of molting
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

### Parameter inputs for probabilit of maturing
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

### Parameter inputs for probabilit of maturing
# EXTRA PARS: Initial  Lower_bound  Upper_bound Prior_type      Prior_1      Prior_2  Phase Reltve 

# Using custom mature probability
#Pre-specified mature probability
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0060 0.0170 0.0510 0.1060 0.1490 0.1780 0.1810 0.1560 0.1930 0.2930 0.4060 0.5320 0.6690 0.8150 0.9140 0.9680 0.9930 0.9870 0.9890 0.9960 
0.0000 0.0000 0.0010 0.0030 0.0220 0.0580 0.0920 0.1240 0.1910 0.2950 0.3550 0.3710 0.3790 0.3770 0.4360 0.5570 0.6930 0.8440 0.9400 0.9800 1.0000 1.0000 
0.0000 0.0000 0.0010 0.0020 0.0170 0.0450 0.0910 0.1540 0.2040 0.2390 0.2680 0.2910 0.3020 0.3000 0.4090 0.6300 0.8000 0.9180 0.9770 0.9760 0.9820 0.9940 
0.0000 0.0000 0.0040 0.0110 0.0580 0.1440 0.2030 0.2360 0.2990 0.3910 0.4520 0.4820 0.5190 0.5650 0.6620 0.8120 0.9130 0.9670 0.9950 0.9980 1.0000 1.0000 
0.0000 0.0000 0.0050 0.0140 0.0370 0.0740 0.1440 0.2470 0.3880 0.5640 0.6630 0.6850 0.7080 0.7320 0.7960 0.9000 0.9600 0.9760 0.9880 0.9960 1.0000 1.0000 
0.0000 0.0000 0.0100 0.0290 0.0580 0.0980 0.1510 0.2170 0.3010 0.4040 0.4770 0.5190 0.5340 0.5220 0.6000 0.7670 0.8830 0.9490 0.9860 0.9950 1.0000 1.0000 
0.0000 0.0000 0.0050 0.0140 0.0410 0.0840 0.1270 0.1680 0.2370 0.3320 0.4240 0.5110 0.5880 0.6550 0.7250 0.7970 0.8710 0.9450 0.9800 0.9740 0.9790 0.9930 
0.0000 0.0000 0.0000 0.0000 0.0350 0.1050 0.1520 0.1750 0.2280 0.3110 0.3730 0.4140 0.4500 0.4810 0.5950 0.7940 0.9080 0.9400 0.9670 0.9890 1.0000 1.0000 
0.0000 0.0000 0.0020 0.0060 0.0340 0.0840 0.1290 0.1700 0.2090 0.2480 0.2650 0.2610 0.2930 0.3620 0.4850 0.6630 0.8070 0.9190 0.9810 0.9940 1.0000 1.0000 
0.0000 0.0000 0.0050 0.0150 0.0460 0.0990 0.1490 0.1940 0.2670 0.3670 0.4080 0.3890 0.4060 0.4600 0.5790 0.7620 0.8820 0.9390 0.9730 0.9860 0.9940 0.9980 
0.0000 0.0000 0.0230 0.0680 0.1170 0.1700 0.2230 0.2760 0.3490 0.4450 0.4690 0.4220 0.4560 0.5720 0.6990 0.8350 0.9170 0.9440 0.9680 0.9890 1.0000 1.0000 
0.0000 0.0000 0.0040 0.0110 0.0280 0.0550 0.1260 0.2420 0.3240 0.3740 0.4010 0.4060 0.4530 0.5450 0.6770 0.8510 0.9460 0.9630 0.9780 0.9930 1.0000 1.0000 
0.0000 0.0000 0.0000 0.0000 0.0330 0.0990 0.1420 0.1610 0.2000 0.2600 0.3320 0.4160 0.4920 0.5590 0.6720 0.8320 0.9230 0.9450 0.9670 0.9890 1.0000 1.0000 
0.0000 0.0000 0.0000 0.0000 0.0070 0.0210 0.0870 0.2030 0.2880 0.3410 0.3590 0.3420 0.3700 0.4430 0.5800 0.7820 0.9130 0.9710 1.0000 1.0000 1.0000 1.0000 
0.0000 0.0000 0.0130 0.0400 0.0720 0.1100 0.1490 0.1890 0.2220 0.2510 0.2910 0.3440 0.3690 0.3660 0.4680 0.6770 0.8200 0.8980 0.9260 0.9040 0.9190 0.9730 
0.0000 0.0000 0.0000 0.0000 0.0370 0.1120 0.1900 0.2700 0.3510 0.4310 0.5330 0.6560 0.6840 0.6170 0.6330 0.7340 0.8280 0.9140 0.9680 0.9890 1.0000 1.0000 
0.0000 0.0000 0.0140 0.0430 0.0780 0.1190 0.1650 0.2160 0.2950 0.4010 0.4850 0.5470 0.5550 0.5080 0.5940 0.8130 0.9330 0.9520 0.9680 0.9810 0.9900 0.9970 
0.0000 0.0000 0.0190 0.0570 0.0950 0.1330 0.2220 0.3600 0.4340 0.4430 0.4420 0.4300 0.4780 0.5880 0.6720 0.7310 0.8080 0.9020 0.9620 0.9870 1.0000 1.0000 
0.0000 0.0000 0.0000 0.0000 0.0320 0.0970 0.2040 0.3550 0.4590 0.5140 0.5310 0.5070 0.4850 0.4650 0.5550 0.7540 0.8680 0.8970 0.8670 0.7790 0.8010 0.9340 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0000 0.0000 0.1000 0.3000 0.3880 0.3630 0.3590 0.3760 0.3480 0.2750 0.2080 0.1490 0.2450 0.4960 0.7160 0.9050 1.0000 1.0000 1.0000 1.0000 
0.0000 0.0000 0.0010 0.0020 0.0250 0.0710 0.1760 0.3420 0.4760 0.5790 0.6240 0.6110 0.6800 0.8330 0.9250 0.9560 0.9640 0.9490 0.9570 0.9860 1.0000 1.0000 
0.0000 0.0000 0.0000 0.0000 0.0310 0.0940 0.1330 0.1460 0.2560 0.4610 0.5440 0.5050 0.4830 0.4810 0.6100 0.8700 1.0000 1.0000 0.9640 0.8930 0.8930 0.9640 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0000 0.0000 0.0160 0.0490 0.1240 0.2430 0.3570 0.4660 0.5300 0.5500 0.6000 0.6800 0.7700 0.8710 0.9320 0.9530 0.9720 0.9910 1.0000 1.0000 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0000 0.0000 0.0310 0.0930 0.1900 0.3210 0.4200 0.4880 0.5340 0.5560 0.5720 0.5810 0.6830 0.8770 0.9810 0.9940 1.0000 1.0000 1.0000 1.0000 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.0090 0.0280 0.1040 0.2370 0.3090 0.3210 0.3550 0.4100 0.4800 0.5630 0.6620 0.7780 0.8650 0.9220 0.9610 0.9800 0.9920 0.9970 1.0000 1.0000 
0.0000 0.0000 0.0080 0.0230 0.0660 0.1360 0.2320 0.3540 0.4420 0.4970 0.5720 0.6650 0.7570 0.8490 0.9180 0.9650 0.9910 0.9970 1.0000 1.0000 1.0000 1.0000 
0.0000 0.0000 0.0420 0.1250 0.1620 0.1540 0.1700 0.2090 0.2440 0.2730 0.3200 0.3860 0.5190 0.7180 0.8350 0.8700 0.9090 0.9520 0.9800 0.9930 1.0000 1.0000 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 
0.0000 0.0000 0.1110 0.3340 0.4770 0.5400 0.5740 0.5790 0.5740 0.5600 0.5610 0.5770 0.6120 0.6640 0.7570 0.8930 0.9710 0.9900 1.0000 1.0000 1.0000 1.0000 
0.0000 0.0000 0.0030 0.0080 0.0370 0.0900 0.1590 0.2450 0.3920 0.5980 0.6630 0.5850 0.5650 0.6030 0.6870 0.8150 0.9000 0.9410 0.9720 0.9910 1.0000 1.0000 
0.0000 0.0050 0.0220 0.0470 0.0810 0.1260 0.1820 0.2480 0.3200 0.3960 0.4710 0.5440 0.6190 0.6950 0.7700 0.8380 0.8960 0.9410 0.9730 0.9940 1.0000 1.0000 
0.0010 0.0150 0.0350 0.0630 0.0990 0.1440 0.1960 0.2560 0.3200 0.3880 0.4590 0.5320 0.6030 0.6700 0.7330 0.7900 0.8430 0.8920 0.9320 0.9630 0.9850 0.9990 
0.0000 0.0000 0.0100 0.0290 0.0670 0.1230 0.1840 0.2500 0.3210 0.3960 0.4450 0.4670 0.5030 0.5510 0.6480 0.7950 0.8930 0.9450 0.9730 0.9770 0.9840 0.9950 

# Using custom mature probability
#Pre-specified mature probability
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 
0.0061 0.0189 0.0584 0.1775 0.4667 0.7761 0.8119 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 


## ==================================================================================== ##
## NATURAL MORTALITY PARAMETER CONTROLS                                                 ##
## ==================================================================================== ##
## 
# Relative: 0 - absolute values; 1+ - based on another M-at-size vector (indexed by ig)
# Type: 0 for standard; 1: Spline
#  For spline: set extra to the number of knots, the parameters are the knots (phase -1) and the log-differences from base M
# Extra: control the number of knots for splines
# Brkpts: number of changes in M by size
# Mirror: Mirror M-at-size over to that for another partition (indexed by ig)
# Block: Block number for time-varying M-at-size
# Block_fn: 0:absolute values; 1:exponential
# Env_L: Environmental link - options are 0: none; 1:additive; 2:multiplicative; 3:exponential
# EnvL_var: Environmental variable
# RW: 0 for no random walk changes; 1 otherwise
# RW_blk: Block number for random walks
# Sigma_RW: Sigma for the random walk parameters
# Mirror_RW: Should time-varying aspects be mirrored (Indexed by ig)
## Relative?   Type   Extra  Brkpts  Mirror   Block  Blk_fn Env_L   EnvL_Vr      RW  RW_blk Sigma_RW Mirr_RW
          0       0       0       0       0       1       1       0       0       0       0   0.3000       0
          1       0       0       0       0       1       1       0       0       0       0   0.3000       0
          0       0       0       0       0       1       1       0       0       0       0   0.3000       0
          3       0       0       0       0       1       1       0       0       0       0   0.3000       0
 # sex*maturity state: male & 1
 # sex*maturity state: male & 2
 # sex*maturity state: female & 1
 # sex*maturity state: female & 2

#      Initial    Lower_bound    Upper_bound  Prior_type        Prior_1        Prior_2  Phase 
    0.27100000     0.15000000     0.70000000           1     0.27100000     0.00454000      4 # M_base_male_mature
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_male_mature_block_group_1_block_1
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_male_mature_block_group_1_block_2
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000     -4 # M_male_mature_block_group_1_block_3
    0.00000000    -1.00000000     1.00000000           0     0.00000000     0.00000000      4 # M_base_male_immature
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000     -4 # M_male_immature_block_group_1_block_1
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_male_immature_block_group_1_block_2
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000     -4 # M_male_immature_block_group_1_block_3
    0.27100000     0.15000000     0.70000000           1     0.27100000     0.00454000      4 # M_base_female_mature
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_female_mature_block_group_1_block_1
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_female_mature_block_group_1_block_2
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000     -4 # M_female_mature_block_group_1_block_3
    0.00000000    -1.00000000     1.00000000           0     0.00000000     0.00000000      4 # M_base_female_immature
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_female_immature_block_group_1_block_1
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000      4 # M_female_immature_block_group_1_block_2
    0.00000000    -1.00000000    10.00000000           0     0.00000000     0.25000000     -4 # M_female_immature_block_group_1_block_3

## ==================================================================================== ##
## SELECTIVITY PARAMETERS CONTROLS                                                      ##
## ==================================================================================== ##
## 
# ## Selectivity parameter controls
# ## Selectivity (and retention) types
# ##  <0: Mirror selectivity
# ##   0: Nonparametric selectivity (one parameter per class)
# ##   1: Nonparametric selectivity (one parameter per class, constant from last specified class)
# ##   2: Logistic selectivity (inflection point and width (i.e. 1/slope))
# ##   3: Logistic selectivity (50% and 95% selection)
# ##   4: Double normal selectivity (3 parameters)
# ##   5: Flat equal to zero (1 parameter; phase must be negative)
# ##   6: Flat equal to one (1 parameter; phase must be negative)
# ##   7: Flat-topped double normal selectivity (4 parameters)
# ##   8: Declining logistic selectivity with initial values (50% and 95% selection plus extra)
# ##   9: Cubic-spline (specified with knots and values at knots)
# ##      Inputs: knots (in length units); values at knots (0-1) - at least one should have phase -1
# ##  10: One parameter logistic selectivity (inflection point and slope)
# ##  11: Pre-specified selectivity (matrix by year and class)
# ##  12: Spline with 0 until one size-class and 1 after another
# ##      Inputs: knots (in length units); values at knots (0-1) - at least one should have phase -1
# ##  13: Stacked logistic
## Selectivity specifications --
# ## Extra (type 1): number of selectivity parameters to estimated
# #  Pot_Fishery Trawl_Bycatch NMFS_Trawl_1982 NMFS_Trawl_1989
 1 0 1 1 # is selectivity sex=specific? (1=Yes; 0=No)
 2 2 0 0 # male selectivity type
 2 2 0 0 # female selectivity type
 0 0 0 0 # selectivity within another gear
 0 0 0 0 # male extra parameters for each pattern
 0 0 0 0 # female extra parameters for each pattern
 1 1 0 0 # male: is maximum selectivity at size forced to equal 1 (1) or not (0)
 1 1 0 0 # female: is maximum selectivity at size forced to equal 1 (1) or not (0)
 0 0 0 0 # size-class at which selectivity is forced to equal 1 (ignored if the previous input is 1)
 0 0 0 0 # size-class at which selectivity is forced to equal 1 (ignored if the previous input is 1)
## Retention specifications --
 1 0 0 0 # is retention sex=specific? (1=Yes; 0=No)
 2 6 6 6 # male retention type
 6 6 6 6 # female retention type
 1 0 0 0 # male retention flag (0 = no, 1 = yes)
 0 0 0 0 # female retention flag (0 = no, 1 = yes)
 0 0 0 0 # male extra parameters for each pattern
 0 0 0 0 # female extra parameters for each pattern
 0 0 0 0 # male - should maximum retention be estimated for males (1=Yes; 0=No)
 0 0 0 0 # female - should maximum retention be estimated for females (1=Yes; 0=No)

## General parameter specificiations 
##  Initial: Initial value for the parameter (must lie between lower and upper)
##  Lower & Upper: Range for the parameter
##  Prior type:
##   0: Uniform   - parameters are the range of the uniform prior
##   1: Normal    - parameters are the mean and sd
##   2: Lognormal - parameters are the mean and sd of the log
##   3: Beta      - parameters are the two beta parameters [see dbeta]
##   4: Gamma     - parameters are the two gamma parameters [see dgamma]
##  Phase: Set equal to a negative number not to estimate
##  Relative: 0: absolute; 1 relative 
##  Block: Block number for time-varying selectivity   
##  Block_fn: 0:absolute values; 1:exponential
##  Env_L: Environmental link - options are 0:none; 1:additive; 2:multiplicative; 3:exponential
##  EnvL_var: Environmental variable
##  RW: 0 for no random walk changes; 1 otherwise
##  RW_blk: Block number for random walks
##  Sigma_RW: Sigma used for the random walk

# Inputs for type*sex*fleet: selectivity male Pot_Fishery
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
           105.711400     5.000000   186.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Pot_Fishery_male_base_Logistic_mean
             4.997241     0.010000    20.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Pot_Fishery_male_base_Logistic_cv

# Inputs for type*sex*fleet: selectivity male Trawl_Bycatch
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
           109.931000     5.000000   185.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Trawl_Bycatch_male_base_Logistic_mean
            11.868260     0.010000    20.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Trawl_Bycatch_male_base_Logistic_cv

# Inputs for type*sex*fleet: selectivity male NMFS_Trawl_1982
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
             0.136100     0.000010     1.000000          1     0.136100     0.117800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_1
             0.248600     0.000010     1.000000          1     0.248600     0.084700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_2
             0.350400     0.000010     1.000000          1     0.350400     0.073800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_3
             0.428600     0.000010     1.000000          1     0.428600     0.072100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_4
             0.474000     0.000010     1.000000          1     0.474000     0.072200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_5
             0.486100     0.000010     1.000000          1     0.486100     0.071200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_6
             0.474000     0.000010     1.000000          1     0.474000     0.069200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_7
             0.452000     0.000010     1.000000          1     0.452000     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_8
             0.432900     0.000010     1.000000          1     0.432900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_9
             0.423100     0.000010     1.000000          1     0.423100     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_10
             0.422300     0.000010     1.000000          1     0.422300     0.067100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_11
             0.427800     0.000010     1.000000          1     0.427800     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_12
             0.438900     0.000010     1.000000          1     0.438900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_13
             0.458600     0.000010     1.000000          1     0.458600     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_14
             0.490600     0.000010     1.000000          1     0.490600     0.067500      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_15
             0.535700     0.000010     1.000000          1     0.535700     0.067700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_16
             0.591800     0.000010     1.000000          1     0.591800     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_17
             0.656900     0.000010     1.000000          1     0.656900     0.068700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_18
             0.732300     0.000010     1.000000          1     0.732300     0.068300      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_19
             0.820900     0.000010     1.000000          1     0.820900     0.070200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_20
             0.923300     0.000010     1.000000          1     0.923300     0.079200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_21
             0.999900     0.000010     1.000000          1     0.999900     0.108700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_male_base_class_22

# Inputs for type*sex*fleet: selectivity male NMFS_Trawl_1989
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
             0.136100     0.000010     1.000000          1     0.136100     0.117800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_1
             0.248600     0.000010     1.000000          1     0.248600     0.084700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_2
             0.350400     0.000010     1.000000          1     0.350400     0.073800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_3
             0.428600     0.000010     1.000000          1     0.428600     0.072100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_4
             0.474000     0.000010     1.000000          1     0.474000     0.072200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_5
             0.486100     0.000010     1.000000          1     0.486100     0.071200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_6
             0.474000     0.000010     1.000000          1     0.474000     0.069200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_7
             0.452000     0.000010     1.000000          1     0.452000     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_8
             0.432900     0.000010     1.000000          1     0.432900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_9
             0.423100     0.000010     1.000000          1     0.423100     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_10
             0.422300     0.000010     1.000000          1     0.422300     0.067100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_11
             0.427800     0.000010     1.000000          1     0.427800     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_12
             0.438900     0.000010     1.000000          1     0.438900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_13
             0.458600     0.000010     1.000000          1     0.458600     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_14
             0.490600     0.000010     1.000000          1     0.490600     0.067500      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_15
             0.535700     0.000010     1.000000          1     0.535700     0.067700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_16
             0.591800     0.000010     1.000000          1     0.591800     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_17
             0.656900     0.000010     1.000000          1     0.656900     0.068700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_18
             0.732300     0.000010     1.000000          1     0.732300     0.068300      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_19
             0.820900     0.000010     1.000000          1     0.820900     0.070200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_20
             0.923300     0.000010     1.000000          1     0.923300     0.079200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_21
             0.999900     0.000010     1.000000          1     0.999900     0.108700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_male_base_class_22

# Inputs for type*sex*fleet: selectivity female Pot_Fishery
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
            74.856720     5.000000   150.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Pot_Fishery_female_base_Logistic_mean
             4.187324     0.010000    20.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Sel_Pot_Fishery_female_base_Logistic_cv

# Inputs for type*sex*fleet: selectivity female NMFS_Trawl_1982
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
             0.136100     0.000010     1.000000          1     0.136100     0.117800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_1
             0.248600     0.000010     1.000000          1     0.248600     0.084700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_2
             0.350400     0.000010     1.000000          1     0.350400     0.073800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_3
             0.428600     0.000010     1.000000          1     0.428600     0.072100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_4
             0.474000     0.000010     1.000000          1     0.474000     0.072200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_5
             0.486100     0.000010     1.000000          1     0.486100     0.071200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_6
             0.474000     0.000010     1.000000          1     0.474000     0.069200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_7
             0.452000     0.000010     1.000000          1     0.452000     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_8
             0.432900     0.000010     1.000000          1     0.432900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_9
             0.423100     0.000010     1.000000          1     0.423100     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_10
             0.422300     0.000010     1.000000          1     0.422300     0.067100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_11
             0.427800     0.000010     1.000000          1     0.427800     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_12
             0.438900     0.000010     1.000000          1     0.438900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_13
             0.458600     0.000010     1.000000          1     0.458600     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_14
             0.490600     0.000010     1.000000          1     0.490600     0.067500      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_15
             0.535700     0.000010     1.000000          1     0.535700     0.067700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_16
             0.591800     0.000010     1.000000          1     0.591800     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_17
             0.656900     0.000010     1.000000          1     0.656900     0.068700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_18
             0.732300     0.000010     1.000000          1     0.732300     0.068300      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_19
             0.820900     0.000010     1.000000          1     0.820900     0.070200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_20
             0.923300     0.000010     1.000000          1     0.923300     0.079200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_21
             0.999900     0.000010     1.000000          1     0.999900     0.108700     -3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1982_female_base_class_22

# Inputs for type*sex*fleet: selectivity female NMFS_Trawl_1989
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
             0.136100     0.000010     1.000000          1     0.136100     0.117800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_1
             0.248600     0.000010     1.000000          1     0.248600     0.084700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_2
             0.350400     0.000010     1.000000          1     0.350400     0.073800      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_3
             0.428600     0.000010     1.000000          1     0.428600     0.072100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_4
             0.474000     0.000010     1.000000          1     0.474000     0.072200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_5
             0.486100     0.000010     1.000000          1     0.486100     0.071200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_6
             0.474000     0.000010     1.000000          1     0.474000     0.069200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_7
             0.452000     0.000010     1.000000          1     0.452000     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_8
             0.432900     0.000010     1.000000          1     0.432900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_9
             0.423100     0.000010     1.000000          1     0.423100     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_10
             0.422300     0.000010     1.000000          1     0.422300     0.067100      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_11
             0.427800     0.000010     1.000000          1     0.427800     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_12
             0.438900     0.000010     1.000000          1     0.438900     0.067200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_13
             0.458600     0.000010     1.000000          1     0.458600     0.067000      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_14
             0.490600     0.000010     1.000000          1     0.490600     0.067500      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_15
             0.535700     0.000010     1.000000          1     0.535700     0.067700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_16
             0.591800     0.000010     1.000000          1     0.591800     0.068200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_17
             0.656900     0.000010     1.000000          1     0.656900     0.068700      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_18
             0.732300     0.000010     1.000000          1     0.732300     0.068300      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_19
             0.820900     0.000010     1.000000          1     0.820900     0.070200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_20
             0.923300     0.000010     1.000000          1     0.923300     0.079200      3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_21
             0.999900     0.000010     1.000000          1     0.999900     0.108700     -3      0      0      0      0      0      0   0.3000 # Sel_NMFS_Trawl_1989_female_base_class_22

# Inputs for type*sex*fleet: retention male Pot_Fishery
# MAIN PARS:  Initial  Lower_bound  Upper_bound Prior_type     Prior_1      Prior_2  Phase  Block Blk_fn  Env_L Env_vr     RW RW_Blk RW_Sigma
            98.039190     1.000000   190.000000          1    98.000000    10.000000      4      0      0      0      0      0      0   0.3000 # Ret_Pot_Fishery_male_base_Logistic_mean
             2.197131     0.001000    20.000000          0     1.000000   999.000000      4      0      0      0      0      0      0   0.3000 # Ret_Pot_Fishery_male_base_Logistic_cv

# pre-specified selectivity/retention (ordered by type, sex, fleet and year)

## ==================================================================================== ##
## CATCHABILITY PARAMETER CONTROLS                                                      ##
## ==================================================================================== ##
## 
# Catchability (specifications)
# Analytic: should q be estimated analytically (1) or not (0)
# Lambda: the weight lambda
# Emphasis: the weighting emphasis
# Block: Block number for time-varying q
# Block_fn: 0:absolute values; 1:exponential
# Env_L: Environmental link - options are 0: none; 1:additive; 2:multiplicative; 3:exponential
# EnvL_var: Environmental variable
# RW: 0 for no random walk changes; 1 otherwise
# RW_blk: Block number for random walks
# Sigma_RW: Sigma for the random walk parameters
## Analytic  Lambda Emphass  Mirror   Block   Env_L EnvL_Vr      RW  RW_blk Sigma_RW
          0       1       1       0       0       0       0       0       0   0.3000 
          0       1       1       0       0       0       0       0       0   0.3000 
          0       1       1       0       0       0       0       0       0   0.3000 
          0       1       1       0       0       0       0       0       0   0.3000 
# Catchability (parameters)
#      Initial    Lower_bound    Upper_bound  Prior_type        Prior_1        Prior_2  Phase 
    1.00000000     0.01000000     1.01000000           0     0.84313600     0.03000000     -5 # Survey_q_parameter_1
    1.00000000     0.01000000     1.01000000           0     0.84313600     0.03000000     -5 # Survey_q_parameter_2
    1.00000000     0.01000000     1.01000000           0     0.45136000     0.50000000     -5 # Survey_q_parameter_3
    1.00000000     0.01000000     1.01000000           0     0.45313600     0.50000000     -5 # Survey_q_parameter_4

## ==================================================================================== ##
## ADDITIONAL CV PARAMETER CONTROLS                                                     ##
## ==================================================================================== ##
## 
# Additiional CV controls (specifications)
# Mirror: should additional variance be mirrored (value > 1) or not (0)?
# Block: Block number for time-varying additional variance
# Block_fn: 0:absolute values; 1:exponential
# Env_L: Environmental link - options are 0: none; 1:additive; 2:multiplicative; 3:exponential
# EnvL_var: Environmental variable
# RW: 0 for no random walk changes; 1 otherwise
# RW_blk: Block number for random walks
# Sigma_RW: Sigma for the random walk parameters
##   Mirror   Block   Env_L EnvL_Vr     RW   RW_blk Sigma_RW
          0       0       0       0       0       0   0.3000 
          0       0       0       0       0       0   0.3000 
          0       0       0       0       0       0   0.3000 
          0       0       0       0       0       0   0.3000 
## Mirror Block Env_L EnvL_Var  RW RW_blk Sigma_RW
# Additional variance (parameters)
#      Initial    Lower_bound    Upper_bound  Prior_type        Prior_1        Prior_2  Phase 
    0.00010000     0.00001000    10.00000000           0     1.00000000   100.00000000     -4 # Add_cv_parameter_1
    0.00010000     0.00001000    10.00000000           0     1.00000000   100.00000000     -4 # Add_cv_parameter_2
    0.00010000     0.00001000    10.00000000           0     1.00000000   100.00000000     -4 # Add_cv_parameter_3
    0.00010000     0.00001000    10.00000000           0     1.00000000   100.00000000     -4 # Add_cv_parameter_4

## ==================================================================================== ##
## CONTROLS ON F                                                                        ##
## ==================================================================================== ##
## 
# Controls on F
#   Initial_male_F  Initial_fem_F   Pen_SD (mal)   Pen_SD (fem) Phz_mean_F_mal Phz_mean_F_fem   Lower_mean_F   Upper_mean_F Low_ann_male_F  Up_ann_male_F    Low_ann_f_F     Up_ann_f_F
          1.000000       0.050500       0.500000      45.500000       1.000000       1.000000     -12.000000       4.000000     -10.000000      10.000000     -10.000000      10.000000  # Pot_Fishery
          0.018000       1.000000       0.500000      45.500000       1.000000      -1.000000     -12.000000       4.000000     -10.000000      10.000000     -10.000000      10.000000  # Trawl_Bycatch
          0.000000       0.000000       2.000000      20.000000      -1.000000      -1.000000     -12.000000       4.000000     -10.000000      10.000000     -10.000000      10.000000  # NMFS_Trawl_1982
          0.000000       0.000000       2.000000      20.000000      -1.000000      -1.000000     -12.000000       4.000000     -10.000000      10.000000     -10.000000      10.000000  # NMFS_Trawl_1989

## ==================================================================================== ##
## SIZE COMPOSITIONS OPTIONS                                                            ##
## ==================================================================================== ##
## 
# Options when fitting size-composition data
## Likelihood types: 
##  1:Multinomial with estimated/fixed sample size
##  2:Robust approximation to multinomial
##  3:logistic normal
##  4:multivariate-t
##  5:Dirichlet
##  6:Dirichlet-Alt (Thorson et al 2016 rec'd)

#  Pot_Fishery Pot_Fishery Pot_Fishery Trawl_Bycatch Trawl_Bycatch NMFS_Trawl_1982 NMFS_Trawl_1989 NMFS_Trawl_1982 NMFS_Trawl_1989 NMFS_Trawl_1982 NMFS_Trawl_1989 NMFS_Trawl_1982 NMFS_Trawl_1989
#   M    M    F    F    M    F    F    M    M    F    F    M    M  
#  ret  tot  dsc  dsc  dsc  dsc  dsc  dsc  dsc  dsc  dsc  dsc  dsc 
#  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S  N+S 
#  I+M  I+M  I+M  I+M  I+M  imm  imm  imm  imm  mat  mat  mat  mat 
      2      2      2      2      2      2      2      2      2      2      2      2      2 # Type of likelihood
      0      0      0      0      0      0      0      0      0      0      0      0      0 # Auto tail compression
 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 # Auto tail compression (pmin)
      1      2      3      4      5      6      7      8      9     10     11     12     13 # Composition aggregator codes
      1      1      1      1      1      2      2      2      2      2      2      2      2 # Set to 1 for catch-based predictions; 2 for survey or total catch predictions
 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 # Lambda for effective sample size
 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 # Lambda for overall likelihood
      0      0      0      0      0      0      0      0      0      0      0      0      0 # Survey to set Q for this comp

# Effective sample size parameters (number matches max(Composition Aggregator code)) 
#      Initial    Lower_bound    Upper_bound  Prior_type        Prior_1        Prior_2  Phase 
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_1(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_2(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_3(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_4(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_5(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_6(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_7(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_8(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_9(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_10(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_11(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_12(possibly extended)
    1.00000000     0.10000000     5.00000000           0     0.00000000   999.00000000     -4 # Overdispersion_parameter_for_size_comp_13(possibly extended)

## ==================================================================================== ##
## EMPHASIS FACTORS                                                                     ##
## ==================================================================================== ##

1.0000 # Emphasis on tagging data

 1.0000 1.0000 1.0000 1.0000 # Emphasis on Catch: (by catch dataframes)

# Weights for penalties 1, 11, and 12
#   Mean_M_fdevs | Mean_F_fdevs |  Ann_M_fdevs |  Ann_F_fdevs
          1.0000         1.0000         0.0000         0.0000 # Pot_Fishery
          1.0000         0.0000         0.0000         0.0000 # Trawl_Bycatch
          0.0000         0.0000         0.0000         0.0000 # NMFS_Trawl_1982
          0.0000         0.0000         0.0000         0.0000 # NMFS_Trawl_1989

## Emphasis Factors (Priors/Penalties: 13 values) ##
 10000.0000	#--Penalty on log_fdev (male+combined; female) to ensure they sum to zero
     0.0000	#--Penalty on mean F by fleet to regularize the solution
     1.0000	#--Not used
     1.0000	#--Not used
    15.0000	#--Not used
     1.0000	#--Smoothness penalty on the recruitment devs
     3.0000	#--Penalty on the difference of the mean_sex_ratio from 0.5
    60.0000	#--Smoothness penalty on molting probability
     3.0000	#--Smoothness penalty on selectivity patterns with class-specific coefficients
     5.0000	#--Smoothness penalty on initial numbers at length
     0.0000	#--Penalty on annual F-devs for males by fleet
     0.0000	#--Penalty on annual F-devs for females by fleet
     0.0000	#--Penalty on deviation parameters

# eof_ctl
9999
