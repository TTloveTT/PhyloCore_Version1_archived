PhyloCore User Manual




 
======================== DEPENDENCY======================== 
PhyloCore depends on Bioperl. Make sure that Bioperl 1.5.2 or later (http://www.bioperl.org/wiki/Getting_BioPerl) is installed.

A script named 'install.pl' is also included with PhyloCore to check and install Bioperl automatically. You need the privilege of the system administrator to run the script. See below for instructions.
 	

======================== INSTALLATION ======================== 

1.Unpack PhyloCore 
	tar -zxvf PhyloCore.tar.gz

2.Run Install (this script will check and install Bioperl, then set the module path in PhyloCore)
	sudo perl install.pl

3.Make PhyloCore executable
	chmod +x PhyloCore.pl


======================== INPUT FILE FORMAT ========================

1.OTU table
      The OTU table should be a tab-delimited table, with sample IDs in the first row, OTU IDs in the first column, and taxonomy in last column. OTU table converted from the biom format can be readily read by PhyloCore (See instructions at http://biom-format.org/).
      
      OTU table format guidelines:
	a. The 1st column header must contain word "OTU".
	b. The last column must contain taxonomy assignment for each OTU. Taxonomic levels should be separated by semicolons. 
		E.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__
	c. OTU table should contain read numbers, but not the relative abundances.


2.OTU tree (optional)
	The OTU tree should be a tree in the newick format. It should contain OTUs in the OTU table.If OTU tree is not provided, PhyloCore will construct a tree with taxa provided in OTU table.	

3.Subset list (optional)
	Specify a sample ID list (with or without group information). Only samples in the list will be used in core identification. This list contain ONLY ONE samples per line, and no header. 



======================== OUTPUT FILES ========================
PhyloCore generates a tab-delimited table, which contains the relative abundance of core taxa in each sample. Core taxa are listed in first column and sample IDs are listed in first row. The last column contains taxonomy for each core taxon.

If user specifies a sample ID list, PhyloCore will also generate a table with only samples in the list.

OTUs appearing in OTU table but not in OTU tree will be written to otu_not_in_tree.log.



======================== RUNNING PhyloCore ========================

Test dataset is provided for troubleshooting.

	RUN test: perl PhyloCore.pl -i test/otu_table.txt [-t test/rep_set.tre] [-s ID_list.txt]

SYNOPSIS
     PhyloCore.pl [options]{-i OTU table}

       [] indicates optional input (order unimportant)

       {} indicates required input (order unimportant)

     Example usage:

	Find Core taxa with >=90% prevalence (default) in all samples
 
 		PhyloCore.pl -i otu_table_taxon.txt 
		
	Find Core taxa with >=90% prevalence (default) in all samples given OTU tree
 
 		PhyloCore.pl -i otu_table_taxon.txt -t otu_rep.tre 

	Find Core taxa with >=80% prevalence in a subset of samples given OTU tree

 		PhyloCore.pl -i otu_table_taxon.txt -t otu_rep.tre -s subsetID_list.txt -p 0.8

       Print help message and exit

       		PhyloCore.pl -h

       Print documentation 
		PhyloCore.pl -man

OPTIONS
       -help or -h
               Print a brief help message and exits.

       -man    Print the manual page and exits.

       
       -t      Specify a OTU tree in the newick format.The OTU tree should contain OTUs in the OTU table.If OTU tree is not provided, PhyloCore will construct a tree with taxa provided in OTU table.

       -s      Specify a sample ID list (with or without group information). Only samples in the list will be used in core identification. The group information, if provided, will be used to identify the weighted core nodes. For weighted core identification, provide a tab-delimited file in which the 1st column is sample ID and the 2nd column is group.

       -o      Specify an output file name. Default: PhyloCore_distribution.txt in current directory

       -prevalence or -p
               Specify a threshold of prevanlence. A node will be considered a core node if its prevalence is above the threshold. Default: 0.9

       -abundance
               Specify a threshold of OTU abundance. OTUs with abundance lower than the threshold in a sample will be considered absent. Default: 0

       -lct_cutoff
              The final taxon for each core will be the lowest common taxonomy shared by the majority of OTUs. The majority threshold is defined by lct_cutoff (0.8 means find the lowest common taxonomy shared by at least 80% of all OTUs in this core node). Default: 0.8

DESCRIPTION
       The program takes a phylogenetic tree and an OTU-table as input and produces a relative abundance table containing all Core taxa.


 
======================== CITATION ======================== 
When publishing work that is based on the results from PhyloCore please cite: ########

======================== LICENSE ========================  
PhyloCore is free software: you may redistribute it and/or modify its under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.

PhyloCore is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).

======================== AUTHOR ========================  
Tiantian Ren, Martin Wu - <http://wolbachia.biology.virginia.edu/WuLab/Home.html>
For any other inquiries send an Email to Tiantian Ren: tr3br@virginia.edu
