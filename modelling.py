"""
High-throughput loop modeling using Modeller. 
Automated loop modelling for non-terminal missing residues given a list of pdb ids.
This tool can also be used for adding missing atoms in PDB structures. 
It is as good as Modeller. Modeller does not give a reliable output when entire domain is 'missing' from the PDB structure.  
"""

import numpy as np
import os
from modeller import *
from modeller.automodel import *
import os
import sys

if not os.path.exists('log_files'):
	os.mkdir('log_files')
if not os.path.exists('fasta'):
	os.mkdir('fasta')
if not os.path.exists('incomplete_structures'):
	os.mkdir('incomplete_structures')


def fasta_reformat(input_seq):
    
    nlines=(len(input_seq)//75)
    rmndr = len(input_seq) % 75
    if rmndr > 0:
        nlines += 1
    
    output_seq = [str(input_seq[i*75:(i+1)*75])+'\n' for i in range(nlines)]
    return output_seq


def write_seq_ali(seqname, new_seq):
    
    fname = '%s.ali' %(seqname)
    f_pml = open(fname, 'a+')
    f_pml.write('>P1;%s \n' %(seqname)) 
    f_pml.write('sequence:%s:::::::0.00: 0.00 \n' %(seqname)) 
    for i in new_seq:
        f_pml.write(i)
    f_pml.close()


def align2d(pid):

	# clean-up PDB file
	cmd='grep -E "ATOM"  incomplete_structures/%s.pdb > %s_incmplt.pdb ' % (pid, pid)

	os.system(cmd)
	env = environ()
	aln = alignment(env)
	mdl = model(env, file='%s_incmplt' %(pid))
	aln.append_model(mdl, align_codes=f'{pid}_incmplt', atom_files=f'{pid}_incmplt.pdb')
	aln.append(file=f'{pid}.ali', align_codes=f'{pid}')
	aln.align2d()
	aln.write(file=f'{pid}_incmplt-{pid}.ali', alignment_format='PIR')

faulty_mdls = []

def loop_modelling(pid):

	# Write new alignment file after removing terminal missing regions
	fname='%s_incmplt-%s-new.ali' %(pid,pid)
	aligned_seq = open(f'{pid}_incmplt-{pid}.ali').readlines()
	pdbseq_end = int(len(aligned_seq)/2)
	refseq_start = int(len(aligned_seq)/2)+3

	pdbseq = ''.join(aligned_seq[3:pdbseq_end]).strip('\n')
	refseq = ''.join(aligned_seq[refseq_start:]).strip('\n')
	
	Nterm = -1
	Cterm = 0

	for i in pdbseq:
		Nterm +=1
		if i !='-':
			break
	
	for i in reversed(pdbseq[:-1]):
		Cterm +=1
		if i !='-':
			break

	pdbseq_new = pdbseq[Nterm:-Cterm]+'*'
	refseq_new = refseq[Nterm:-Cterm]+'*'
	pdbseq_new = pdbseq_new.replace('\n', '') 
	refseq_new = refseq_new.replace('\n', '')

	# generate a list of pdb structures with long missing loops or those with an entire domain missing from the structure. These structures require manual remodeling of missing loops or use of Alpha-fold for building reliable structural models.

	counter = pdbseq_new.count('-')
	if counter>20:
		faulty_mdls.append(pid)
		print(pid)
		
	f_pml = open(fname, 'a+')

	for i in aligned_seq[:3]:
		f_pml.write(i)
	for i in fasta_reformat(pdbseq_new):
		f_pml.write(i)   
	for i in aligned_seq[int(len(aligned_seq)/2):int(len(aligned_seq)/2)+3]:
		f_pml.write(i)
	for i in fasta_reformat(refseq_new):
		f_pml.write(i)
	f_pml.close()

	env = environ()
	a = automodel(env, alnfile=fname,
              knowns=f'{pid}_incmplt', sequence=f'{pid}',
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341)) 
	a.starting_model = 1
	a.ending_model = 1
	a.make()
	
	os.rename(f"{pid}.B99990001.pdb",f"{pid}.pdb")
	cmd='mv %s*0001 log_files' %(pid)
	os.system(cmd)

if __name__ == "__main__":

	try:
		pdbidlist = sys.argv[1]
	except IndexError:
		print("Provide a file with list of pdb-ids (in lower-case)")
		print("Usage: <Modeler-installation-path>/bin/modpy.sh  python modelling.py <pdb-id-list file>")
		sys.exit(0)
	
	pdbids = np.genfromtxt(pdbidlist, dtype=str)

	for p in pdbids:

		# Download PDB structure
		cmd = "wget -O incomplete_structures/%s.pdb  https://files.rcsb.org/download/%s.pdb" %(p, p)
		os.system(cmd)
	
		# Download & read fasta sequence from RCSB PDB
		cmd = "wget -O fasta/%s.fasta https://www.rcsb.org/fasta/entry/%s" %(p, p)
		os.system(cmd)
		input_seq = open(f'fasta/{p}.fasta', 'r').readlines()
		#input_seq = input_seq[1:][0].replace('\n','*')
		input_seq = '/'.join([input_seq[2*n+1:2*n+2][0] for n in range(len(input_seq)//2)]).replace('\n','') + '*'

		# Reformat fasta sequence to Modeler-compatible format & save it as .ali sequence file
		new_seq = fasta_reformat(input_seq)
		write_seq_ali(p, new_seq)

		# Align refromatted sequence with the (incomplete) sequence derived from the PDB structure 
		align2d(p)

		# Remove terminal missing regions and model non-terminal missing loops	
		loop_modelling(p)


np.savetxt('faulty_models.txt', np.array(faulty_mdls), fmt="%s")	
