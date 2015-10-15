#!/usr/bin/perl -w
#use lib "bioperl_home_dir";
use strict;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Tree;
use Bio::Tree::Node;
use Getopt::Long;
use Pod::Usage;


######## Usage #####################
my ($otutable,$otutree,$samplelist,$help,$man) = undef;
my $prevalence_cutoff = 0.9;
my $abundance_cutoff = 0;
my $LCA_cutoff = 0.9;
my $output1 = "PhyloCore_distribution.txt";


pod2usage(-verbose => 1) && exit unless @ARGV;

GetOptions(
    'i=s'=> \$otutable,
    't=s'=>\$otutree,
    's=s'=>\$samplelist,
    'o=s'=>\$output1,
    'prevalence|p=f'=>\$prevalence_cutoff,
    'abundance=f'=>\$abundance_cutoff,
    'lct_cutoff=f'=>\$LCA_cutoff,
    'help|h' => \$help,
    'man'  => \$man,

)  or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $help;
pod2usage(-verbose => 2) && exit if defined $man;
pod2usage(-verbose => 1) && exit unless $otutable;


######### Declare useful variables ######################
#table: otu table hash of hashes ; Taxon:otu->taxon pair; core_nodes: core node internalID ->core node pair  ; core_taxa: hash of hashes, {core taxon}{OTU}= taxon ; samples: sampleIDs ;groups: {group}{sampleID} = 1

my (%table,%Taxon,%reads,%core_nodes,%core_taxa,%samples,%groups);  


#######read in OTU table and save in hash of hashes(%table).######  
open (IN, $otutable) || die "Can't open $otutable\n ";

my $header = undef;
my $flag = 0;
while (<IN>) {
    chomp;
    if (/OTU/){
	$header = $_;
	$flag = 1;
	last;
    }
}
unless ($flag){
    die "OTU table does not have correct header! (e.g OTUID, sampleID1, sampleID2 ..., sampleIDn,taxonomy)";
}

my @sampleID = split /\t/, $header;
pop @sampleID;
shift @sampleID;

while (<IN>){
    chomp;
    my @data=split /\t/, $_;
    my $OTUID = shift @data;
    my $taxon = pop @data;
    unless ($taxon =~ /[A-Za-z]/){
	die "Taxonomy column needed for OTU $OTUID!";
    }
    for my $i (0..$#sampleID){
	if ($data[$i]=~ /^\d+\.0$/ or $data[$i]=~ /^\d+$/ ){
	    $table{$OTUID}{$sampleID[$i]}=$data[$i];	
	    $Taxon{$OTUID}=$taxon;
	    $reads{$sampleID[$i]}+=$data[$i];
	    $samples{$sampleID[$i]}=1;
	    $groups{"all"}{$sampleID[$i]}=1;
	}
	else{
	    die "OTU table can only contain non-negative integer, $OTUID contain illegal number!";
	}
    }
}
close IN || die $!;

########read in a list of sample IDs user wants to test.########
if (defined $samplelist){
    open (IN, $samplelist) || die "Can't open $samplelist\n";
    my (%subset,%subgroups) =();
    while(<IN>){
	chomp;
	my ($s,$g) = split/\t/,$_;
	if (exists $samples{$s}){
	    if (defined $g){	
		$subgroups{$g}{$s} = 1;
	    }
	    else{
		$subgroups{"all"}{$s} = 1;
	    }
	    $subset{$s} = 1;
	}
	else{
	    die "$s is not a correct sample ID!";
	}
    }
    %samples = %subset;
    %groups = %subgroups;
    close IN || die $!;
}


######################## Traverse the tree and find core nodes ##############################
if (!$otutree){

######## if OTU tree is not provided, generate tree based on taxonomy information in OTU table##########

    #Make taxa tree
    my $tree = make_taxa_tree(%Taxon);
    
    #Loop through the tree and find cores based on prevalence cutoff, this step will save the final core nodes in the hash %core_nodes
    my $rootnode =$tree->get_root_node;
    traverse($rootnode); 
    foreach my $nodeinterid (keys %core_nodes){
	#find the full lineage for each core node
	my @lineage = $tree->get_lineage_nodes($core_nodes{$nodeinterid});
	my $full_taxa = undef;
	for my $l (@lineage){
	    $full_taxa .= $l->id.";";      
	}
	$full_taxa .= $core_nodes{$nodeinterid}->id;
	
	#save core info in core_taxa: hash of hashes, {core taxon}{OTU}= taxon   
	if (!$core_nodes{$nodeinterid}->is_Leaf){  #if core is not OTU level
	    foreach my $offspring ($core_nodes{$nodeinterid}->get_all_Descendents()){ #get all leaves, namely OTU ID
		if ($offspring->is_Leaf){
		    $core_taxa{$core_nodes{$nodeinterid}->id}{$offspring->id} = $full_taxa;
		}
	    }
	    
	}
	else {   #if core is OTU level
	    $core_taxa{$core_nodes{$nodeinterid}->id}{$core_nodes{$nodeinterid}->id} = $Taxon{$core_nodes{$nodeinterid}->id};
	}
    }
}


else{

######## if OTU tree is provided, use OTU tree ###########

    #read in OTU newick tree and find root node
    my $tree_input = new Bio::TreeIO(-file => $otutree, -format => "newick"); 
    
    my $tree = $tree_input->next_tree;
    
    unless ($tree) { die "Tree file is not a newick tree!\n";}
    
    my $rootnode =$tree->get_root_node;
    
    #validation: OTUs in OTU table also need to be in the tree, if not, OTU ID will be write to otu_not_in_tree.log
    my %treeOTU =();
    open (ERR, ">otu_not_in_tree.log") || die "Can't otu_not_in_tree.log\n";
    for my $offspring ($rootnode->get_all_Descendents()){
	if ($offspring ->is_Leaf){
	    $treeOTU{$offspring->id} =1;
	}
    }
    for my $otu (keys %table){
	unless (exists $treeOTU{$otu}){
	    print ERR "$otu is not in $otutree!\n";
	}
    }
    close ERR ||die $!;
    


    #loop through the tree and find cores based on prevalence cutoff, this step will save the final core nodes in the hash %core_nodes
    traverse($rootnode);
    
    #find all OTUs within one core, find last common ancestor of all OTUs belong to the core node (except OTU level core), save all core taxa in %core_taxa.
    foreach my $nodeinterid (keys %core_nodes){ 
	my %core_OTUs = ();
	if (!$core_nodes{$nodeinterid}->is_Leaf){  #if core is not OTU level, then find the Last common ancestor.
	    foreach my $offspring ($core_nodes{$nodeinterid}->get_all_Descendents()){ #get all leaves, namely OTU ID
		if ($offspring->is_Leaf && exists $Taxon{$offspring->id}){
		    $core_OTUs{$offspring->id} = $Taxon{$offspring->id};
		}
	    }
	    my $taxa_tree = make_taxa_tree(%core_OTUs);
	    findLCA($taxa_tree);
	}
	else {   #if core is OTU level
	    $core_taxa{$core_nodes{$nodeinterid}->id}{$core_nodes{$nodeinterid}->id} = $Taxon{$core_nodes{$nodeinterid}->id};
	}
    }
}



################### output  ############################

#output OTU distribution in all samples
open (OUT, ">$output1") ||die "Can't open $output1\n";
my $output_all = undef;
$output_all.="Cores\t".join("\t",@sampleID)."\tTaxonomy\t"."\n";
for my $taxon (sort keys %core_taxa){
    next if ($taxon =~ /Bacteria/);
    $output_all.=$taxon."\t";
    for my $i (0..$#sampleID){
	my $sum = 0;
	for my $OTU (keys %{$core_taxa{$taxon}}){
	    $sum += $table{$OTU}{$sampleID[$i]};
	}
	$output_all.=$sum/$reads{$sampleID[$i]}."\t";
    }
    my (undef,$value)= %{$core_taxa{$taxon}};
    $output_all.=$value."\n";    
}
print OUT $output_all;
close OUT ||die $!;


#output OTU distribution in subset, if list provided.
if (defined $samplelist){
    my $output2 = "subset_".$output1;
    open (OUT, ">$output2") ||die "Can't open $output2\n";
    my $output_subset = undef;
    $output_subset.="Cores\t";
    for my $s (keys %samples){
	$output_subset.= $s."\t";
    }
    $output_subset.= "Taxonomy\n";
    for my $taxon (sort keys %core_taxa){
	next if ($taxon =~ /Bacteria/);
	$output_subset.= $taxon."\t";
	for my $s (keys %samples){
	    my $sum = 0;
	    for my $OTU (keys %{$core_taxa{$taxon}}){
		$sum += $table{$OTU}{$s};
	    }
	    $output_subset.= $sum/$reads{$s}."\t";
	}
	my (undef,$value)= %{$core_taxa{$taxon}};
	$output_subset.= $value."\n";    
    } 
    print OUT $output_subset;
    close OUT ||die $!;
}
print "Done!\n";
    

############# subroutines #####################


#subroutine1 calculate prevalence(%) of an array of OTU IDs,given list of subset and abundance cutoff 
sub prevalence {
    my @testOTU= @_;
    my (%exist, $preval)=();
    for my $group (keys %groups){
	$exist{$group} = ();
	for my $OTU (@testOTU) { #ignore OTUs in the tree but not in OTU table
	    if (exists $table{$OTU}){
		#Weight samples by group, prevalance = [exist(group A)/n(group A) + ...+ exist(group Z)/n(group Z)]/ # of groups
		for my $sample (keys %{$groups{$group}}){ #only consider samples in subset list, calulate by group
		    if ($table{$OTU}{$sample}/$reads{$sample}>$abundance_cutoff){   #only consider abundant OTUs
			$exist{$group}{$sample} = 1;
		    }
		}	
	    }
	}
    }
    for my $group (keys %exist){
	$preval += (scalar keys %{$exist{$group}})/(scalar keys %{$groups{$group}});
    }        
    return $preval/(scalar keys %groups); #calculate prevalence  
}

#subroutine2 traverse the tree, recursion, find core
sub traverse {
    my ($node) = @_;
    my @leaves =(); 
     if (!$node->is_Leaf){
	foreach my $offspring ($node->get_all_Descendents()){ #get all leaves, namely OTU ID
	    if ($offspring->is_Leaf){
		push @leaves, $offspring->id;
	    }
	}
    }
    else{
	push @leaves,$node->id;
    }
    if (prevalence(@leaves)>=$prevalence_cutoff){ 	#test condition here, prevalence cutoff
	$core_nodes{$node->internal_id}=$node;
	if ($node->ancestor){
	    delete $core_nodes{($node->ancestor)->internal_id};
	}		
	foreach my $child ($node->each_Descendent){					
	    traverse($child);
	}					
    }	
}

#subroutine3 make tree from taxonomy info of a list of OTUs
sub make_taxa_tree{
    my %list = @_;
    #make tree from taxonomy info of a list of OTUs
    my $tree = Bio::Tree::Tree->new(-root => Bio::Tree::Node->new(-id => "k__Bacteria", -description => 0));  #use description to save how many leaves in each node

    for my $otu (keys %list){
	my @levels = split(/;\s*/,$list{$otu});   #delimiter is (;+ any space) or ;
	my $current_node = $tree->get_root_node();
	for my $i (1 .. $#levels){
	    if ($levels[$i] =~ /__($|[a-z]$)|uncultured|norank/i){ #undefined taxa can be g__ in greengene, __g... and __uncultured in SILVA, ignore those.
		last;
	    }
	    my @offsprings = $current_node->each_Descendent;
	    my $exist = 0;
	    for my $offspring (@offsprings){
		if ($offspring->id eq $levels[$i]){
		    $exist = 1;
		    $current_node->description($current_node->description()+1);
		    $current_node = $offspring;
		    last;
		}
	    }
	    unless ($exist) {
		my $new_node = Bio::Tree::Node->new(-id => $levels[$i],-description => 0);
		$current_node->add_Descendent($new_node);
		$current_node->description($current_node->description()+1);
		$current_node = $new_node;				
	    }
	}		
	my $new_node = Bio::Tree::Node->new(-id => $otu, -description => 0);
	$current_node->add_Descendent($new_node);
	$current_node->description($current_node->description()+1);
	
    }
    return $tree

}

#subroutine4 find lowest common taxonomy of a list of OTUs with taxa tree
sub findLCA {
    my ($tree) = @_;
    
    #loop tree and find lowest common taxonomy
    my $current_node = $tree->get_root_node();
    my $totalOTU = $current_node->description();
    while(!$current_node->is_Leaf){
	my @offsprings = $current_node->each_Descendent;
	my $flag = 0;
	for my $offspring (@offsprings){
	    if ($offspring->description()>= $LCA_cutoff*$totalOTU){
		$current_node = $offspring;
		$flag = 1;
		last;
	    }		
	}
	unless ($flag){   
	    last;	
	}
    }
    my @lineage = $tree->get_lineage_nodes($current_node);
    my $taxon = undef;
    for my $l (@lineage){
	$taxon .= $l->id.";";      
    }
    $taxon .= $current_node->id;
#    if (exists $core_taxa{$current_node->id}){
#	$current_node->id($current_node->id."+");
#    }
    for my $offspring ($current_node->get_all_Descendents){
	if ($offspring->is_Leaf){
	    $core_taxa{$current_node->id}{$offspring->id} = $taxon; #core_taxa: hash of hashes, {core taxon}{OTU}= taxon 
	}
    }
}



    

__END__

=head1 NAME

PhyloCore.pl - Find the core taxa in communities


=head1 SYNOPSIS

PhyloCore.pl [options]{-B<i> OTU table}

{} indicates required input (order unimportant)

[] indicates optional input (order unimportant)

=====File format=====

1.OTU table

    The OTU table should be a tab-delimited table, with sample IDs in the first row, OTU IDs in the first column, and taxonomy in last column. OTU table converted from the biom format can be readily read by PhyloCore (See instructions at http://biom-format.org/). You can find an example in the test directory. 
      
      The followings details the format guidelines:

	a. The 1st column header must contain word "OTU".

	b. The last column must contain taxonomy assignment for each OTU. Taxonomic levels should be separated by semicolons.
 
		E.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__

	c. OTU table should contain read numbers, but not the relative abundances.

2.OTU tree (optional)

	The OTU tree should be a tree in the newick format. It should contain OTUs in the OTU table.If OTU tree is not provided, PhyloCore will construct a tree with taxa provided in OTU table.

3.subset list (optional)

	Users can provide a sample ID list (with or without group information). 
        Only samples in the list will be used in core identification. 
        This list contain ONLY ONE samples per line, and no header. 

        The group information, if provided, will be used to identify the weighted core nodes.
        For weighted core, provide a tab-delimited file in which the 1st column is sample ID and the 2nd column is group.



=====Example usage=====

Print help message and exit

    PhyloCore.pl -h

Find Core taxa with >=90% prevalence (default) in all samples 
 
    PhyloCore.pl -i otu_table_taxon.txt 

Find Core taxa with >=90% prevalence (default) in all samples given OTU tree
 
    PhyloCore.pl -i otu_table_taxon.txt -t otu_rep.tre 

Find Core taxa with >=80% prevalence in a subset of samples given OTU tree

    PhyloCore.pl -i otu_table_taxon.txt -t otu_rep.tre -s subsetID_list.txt -p 0.8

Print documentation 
    PhyloCore.pl -man

=head1 OPTIONS

=over 8

=item B<-help> or B<-h>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-t>

Specify a OTU tree in the newick format.The OTU tree should contain OTUs in the OTU table.If OTU tree is not provided, PhyloCore will construct a tree with taxa provided in OTU table.

=item B<-s>

Specify a sample ID list (with or without group information). Only samples in the list will be used in core identification. The group information, if provided, will be used to identify the weighted core nodes. For weighted core, provide a tab-delimited file in which the 1st column is sample ID and the 2nd column is group.

=item B<-o>

Specify an output file name. Default: PhyloCore_distribution.txt in current directory

=item B<-prevalence> or B<-p>

Specify a threshold of prevanlence. A node will be considered as a core node if its prevalence is above the threshold. Default: 0.9

=item B<-abundance>

Specify a threshold of OTU abundance. OTUs with abundance lower than the threshold in a sample will be considered absent. Default: 0

=item B<-lct_cutoff>

The final taxon for each core will be assigned based on the taxonomy of major (defined by lct_cutoff, 0.8 means 80%) OTUs' taxa in this core. Default: 0.8

=back

=head1 DESCRIPTION  

PhyloCore takes a phylogenetic tree and an OTU-table as input and produces a relative abundance table containing all Core taxa. PhyloCore generates a tab-delimited table, which contains the relative abundance of core taxa in each sample. Core taxa are listed in first column and sample IDs are listed in first row. The last column contains taxonomy for each core taxon. 
If user specifies a sample ID list, PhyloCore will also generate a table with only samples in the list. 


=head1 LICENSE

PhyloCore is free software: you may redistribute it and/or modify its under the terms of the 
GNU General Public License as published by the Free Software Foundation; either version 2 of
the License, or any later version.

PhyloCore is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
GNU General Public License for more details. See L<http://www.gnu.org/licenses/>.

=head1 AUTHOR

Tiantian Ren (tr3br@virginia.edu), Martin Wu - L<http://wolbachia.biology.virginia.edu/WuLab/Home.html>

When publishing work that is based on the results from PhyloCore please cite:

=cut
