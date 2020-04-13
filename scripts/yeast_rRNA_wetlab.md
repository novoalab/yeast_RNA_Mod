# Yeast Ribosomal RNA 

## Producing the data:

### Libraries
For the Weizmann samples, we received total RNA isolated by Trizol extraction\
We further processed the sample by Turbo DNAse treatment with consequent cleanup\
Furthermore, we added PolyA tail to the RNA samples using E.Coli Poly(A) Polymerase


For the sucrose gradient samples, we extracted different fractions and isolated only RNA using Trizol

A classic protocol for the processing of the samples are as follows

#### (Optional) Thawing the total RNA samples and DNAse treatment (Turbo DNAse)
	- This is for Weizmann samples, where we now that they are not treated with DNAse
  * Thaw the RNA samples on ice (Stock concentration around 6 ug/ul)\
  * Once thawed, make aliquots of each tube into 4-5 tubes to avoid freeze-thawing all the RNA\
  * Also make 1:10 dilution into another tube to use for Qubit,TapeStation etc. (Keep this one on ice)\
  * Take 1ul (around 6 ug) and use for Turbo DNAse


In a PCR tube
  * (Manufacturer recommends to dilute nucleic acid concentration to 200 ng/ul, that’s why we are DNAse treating in 50 ul volume)
  * Normal recommendation of DNAse amount is 1 ul, but we do a vigorous treatment to avoid longer incubation

| Reagent	|Initial Conc|Volume|	Final Conc	|
|----------|--------|-----|------|
| Total RNA	||	1 ul	| Max 200 ng/ul|	
| Turbo DNAse||		2 ul 	|	
| Rnase Inhibitor	||	1 ul	|	
| Turbo DNAse Buffer	|10X |	5 ul|	1 X|
| Water		||41 ul|		
| Total		||50 ul	|	

  * Incubate 10 minutes at 37 C 

#### Cleanup using RNA Clean XP Beads

  * Mix the samples with the appropriate volume of beads (1.8X keeps everything above 100bp)
  * Mix the beads with sample 10 times by pipetting (You could also put it on Hula Mixer)
  * Incubate 5-10 minutes at room temperature
  * Spin down the tube and place it on the magnet
  * Remove the supernatant
  * Add 70% freshly prepared ethanol to the tube (200 ul)
  * Incubate for 30 seconds at room temperature
  * Remove the ethanol completely by spinning down and placing back it on magnet.
  * Repeat the washing if you want to get rid of non-specific things more
  * Air-dry the pellet up to 2 minutes (Otherwise they dry out and crack!!!) 
  * Resuspend the beads in water
  * Incubate  5-10 minutes in RT
  * Place the beads on magnet
  * Remove the elute and place it into another tube



#### **(Common steps starting from here)** PolyA tailing using E. Coli Poly(A) Polymerase


|Component	|Volume|
|-------|------|
|RNA + Water|	1ug|
|10X buffer	|2ul|
|ATP(10mM)	|2ul|
|E.coli PAP	|2ul*(Normally 1ul)|
|Rnase Inhibitor|	0.5 ul|
|Total	|20 ul|

  * Incubate 15 minutes at 37 C

#### Cleanup using Zymo RNA Clean and Concentrator


  * Add 2 volumes of RNA Binding Buffer to each sample1 and mix. Example: Mix 100 μl buffer and 50 μl sample.	
  * Add an equal volume of ethanol (95-100%) and mix. Example: Add 150 μl ethanol.
  * Transfer the sample to the Zymo-SpinTM IC Column in a Collection Tube and centrifuge for 30 seconds. Discard the flow-through.
  * Add 400 μl RNA Prep Buffer to the column and centrifuge for 30 seconds. Discard the flow-through.
  * Add 700 μl RNA Wash Buffer to the column and centrifuge for 30 seconds. Discard the flow-through.
  * Add 400 μl RNA Wash Buffer to the column and centrifuge for 2 minutes to ensure complete removal of the wash buffer. Transfer the column carefully into an RNase- free tube (not provided).
  * Add 15 μl DNase/RNase-Free Water directly to the column matrix and centrifuge for 30 seconds.
  * Alternatively, for highly concentrated RNA use ≥6 μl elution. 
  * Also, for even more concentrated RNA, after adding 6 ul water, add 6 ul more and pool those 12 ul. Or if you are aiming for a more more concentrated RNA, take the first elute and reapply.



#### Direct RNA Sequencing Library Prep with Barcoding

  * We will use the barcodes in order to multiplex the libraries

##### Pre-annealing of Oligos

  * Pre-anneal the oligo pairs before the ligation reaction


|Reagent|	Initial Concentration|	Volume (uL)|	Final Conc|
|-------|-------------|------------|----------|
|Oligo A	|100 uM	|1.05|	1.4 uM |
|Oligo B 	|100 uM	|1.05	|1.4 uM|
|Tris pH7.5	|0.1 M|	7.5	|0.01 M|
|NaCl	|0.5 M|	7.5|	0.05 M|
|RNAseOut	||	1|	|
|dH2O	||	50.9|	|
|Total	||	75|	|

 * Heat the mixture for 94°C for 5 mins and ramp down to RT at 0.1°C/s (in PCR machine).


##### Library Prep

  * For each reaction, we will scale down it by half 
  * We will use 250 ng input for unfragmented (around 2000nt average length)


  * In a 0.2 mL thin-walled PCR tube, mix the reagents in the following order (MASTER MIX)
	 
	|Reagent	|Volume (uL)|	8.5 X|
	|-------| -------|-------|
	|NEBNext Quick Ligation Reaction Buffer|	1.5|	12.75|
	|T4 DNA Ligase|	0.75|	6.375|
	|RNA+Water MIxes	|4.5	||

  * Distribute this master mix to each tube
  * Add pre-annealed oligos to each tube (0.5 ul)
  * Mix by pipetting and spin down.
  * Incubate the reaction for 10 mins at RT. 
  * Mix the following reagents together to make the reverse transcription master mix:



|Reagent	|Volume (uL)	|
|-------|-------|
|Nuclease-free water|6.5|
|10mM dNTPs|	1	|
|5x Maxima RT Buffer|	4|	
|Total	12.5	|

  * Add the master mix to the 0.2 mL PCR tube containing the RT adaptor ligated RNA from the “RT Adapter ligation” step. Mix by pipetting.
  * Add 1 uL Maxima RT reverse transcriptase to the reaction and mix by pipetting. 
  * Place the tube in a thermal cycler and incubate at 60°C for 30 mins, and bring the sample to 4°C before proceeding to the next step.  (No heat inactivation)
  * Transfer the sample to a 1.5 mL DNA LoBind Eppendorf tube.
  * Resuspend the stock of Agencourt RNAClean XP beads by vortexing. 
  * Add 36 ul (1.8X) beads to each tube and  mix by pipetting. 
  * Incubate on Hula mixer for 5 mins at RT. 
  * Prepare 2mL of fresh 70% ethanol in nuclease free water. 
  * Spin down the sample and pellet on a magnet. 
  * Keep the tube on the magnet and wash the beads with 200 uL  70% ethanol without disturbing the pellet. (No resuspension on washing buffer). 
  * Just incubate 30 seconds on magnet. 
  * Remove the 70% ethanol and discard. Spin down tubes place back on magnet and remove residual ethanol. 
  * Remove the tube from the magnetic rack and resuspend pellet in 5 uL nuclease-free water. Incubate for 5 mins at RT. 
  * Pellet beads on magnet until the eluate is clear and colourless. 
  * Pipette 5 uL of eluate into a clean 1.5 mL Eppendorf DNA LoBind tube. 
  * Measure cDNA and RNA on Qubit high sensitivity. 
  * RNAs are usually too low
  * Calculate based on cDNA!!!

 Final amount for the library : 200 ng, so 50 ng each



Second Ligation Step

  * In a clean 1.5 mL Eppendorf DNA LoBind tube, mix the reagents in the following order

|Reagent	| Volume (uL)|
|-------|------|
|Reverse-transcribed RNA from the “Reverse Transcription” step|	20|
|NEBNext Quick Ligation Reaction Buffer	|8|
|RNA Adaptor (RMX)	|6|
|Nuclease-free water|	3|
|T4 DNA ligase|	3|
|Total|	40|

  * Mix by pipetting. Incubate for 10 minutes at RT. 
  * Re-suspend the stock of Agencourt RNA Clean XP beads by vortexing, 
  * Add 40 uL of beads to the adaptor ligation reaction and mix by pipetting.
  * Incubate on a Hula mixer for 5 mins at RT. 
  * Spin down the sample and pellet on a magnet. Keep the tube on a magnet and pipette of the supernatant. 
  * Add 150 uL of the Wash Buffer (WSB) to the beads. Resuspend the beads by flicking the tube. Return the tube on the magnetic rack, allow beads to pellet and pipette of the supernatant. Repeat. Allow to air dry for 2 mins.
  * Remove the tube from the magnetic rack and resuspend the pellet in 21 uL Elution Buffer. Incubate for 10 mins at RT. 
  * Pellet beads on magnet until eluate is clear and colourless. 
  * Remove and retain 21 uL of eluate into a clear 1.5 mL Eppendorf DNA LoBind tube. 
  * Measure cDNA on Qubit high sensitivity. 
		○ Library 5 ng/ul 

  * Mix the libraries with 17.5 ul water
  * Add 37.5 ul RRB Buffer and mix well


Flow Cell Priming
	
  * Take out the flowcell from the fridge
  * Insert it to the MinION
  * QC the flowcells :
	- Open MinKNOW software 
	- Select FlowCell type (106)  
	- Click on the flow cell  
	- Click "Check Flow Cell"  
	- Start. 
	

  * Open the priming port by sliding it down clock-wise
  * Adjust 1000ul pipette to 200 ul
  * Put the tip inside the priming port
  * While the tip is in the port, increase the volume of pipette from 200ul up to 220ul 
  * Remove the tip from the port once you see a yellow liquid
  * Mix 600 ul RRB (Mixed and vortexed) with 600 ul water
  * Take 800 ul RRB+water mix and introduce it to priming port:
	- To avoid bubbles, when the tip is really close to entering the priming port, put the first drop outside and then insert the tip to the port. 
	- Insert all the mix slowly to the port AND put the last drop also outside
	- Incubate 5 minutes while priming port is still open
	- Open the SAMPLE PORT 
	- Add 200 ul more RRB+water mix the same way with 1000 ul pipette to the PRIMING PORT AND SEE THE BUBBLES COMING OUT OFF OF THE SAMPLE PORT
	- Now put ONLY ONE DROP of the library through the SAMPLE PORT and start the sequencing:
  * On the flow cell
		- Close both priming and sample ports
		- Close the lid
  * On the software: 
		- New experiment
		- Paste the Run ID 
		- Turn off Live Base Calling
		- Select RNA002
		- Start sequencing
  * Once you are sure that flow cell is good, apply the rest of the library by opening both priming and sample port and loading the library through the sample port. Close the ports again.
