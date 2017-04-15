from django.shortcuts import render, HttpResponse
import string, random
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
import pylab
from primer3 import bindings

tree_map={}

#ACCESSORIES FUNCTIONS to support the views

'''
#RETURNS the Nucleotide Base which is seen the highest number of times in each
in the alignment
##Also substitutes '-' in the alignment with 'N' to prepare for Primer Design
'''
def countIndBase(col, numRows):
    count = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0, 'N' : 0, '-':0, 'W':0,'S':0,'H':0, 'Y':0,'M':0,'R':0}
    for letter in col:
        count[letter] += 1


    baseNuc= max(count, key=count.get)
    if (baseNuc=='-'):
        return 'N'
    #print " The returning base is: ", baseNuc
    else:
        return max(count, key=count.get)

#Takes FILENAME as INPUT and performs MSA using CLUSTALW through COMMANDLINE
def multipleSeqAlignment(fileName):

    cline=ClustalwCommandline("clustalw", infile=fileName, output="fasta", seed=100)
    cline()


def node_name(leaf):
    if leaf.name !=None:
        return tree_map[leaf.name]
    else:
        return leaf.name

def tree_draw_newick():

    multipleSeqAlignment("input.fasta")

    alignment= AlignIO.read("input.fasta", "fasta")
    number_col=alignment.get_alignment_length()
    number_row=len(alignment)


    consensusSeq=''
    for each in list(range(number_col)):
        column=alignment[:,each]

        consensusSeq=consensusSeq + countIndBase(column, number_row)
    #return HttpResponse(consensusSeq)
    seq_args={'SEQUENCE_ID':'CONSENSUS_SEQ', 'SEQUENCE_TEMPLATE': consensusSeq}
    global_args= {'PRIMER_OPT_SIZE': 20,'PRIMER_PICK_INTERNAL_OLIGO': 1,'PRIMER_PRODUCT_SIZE_RANGE': [[175,200],[200,225]],}


    results=bindings.designPrimers(seq_args, global_args)

    treeNewick= Phylo.read("input.dnd", "newick")
    #Phylo.draw_ascii(treeNewick)
    treeNewick.rooted=True
    #note sure what this does
    treeNewick.ladderize()
    Phylo.draw(treeNewick, do_show=False, label_func= node_name)
    pylab.axis('on')
    pylab.savefig("/home/savill88/Desktop/senior-project/PriMSA/django_PriMSA/primsa/static/images/phylotree", format='png', bbox_inches='tight', dpi=250)




    return results

#END of ACCESSORIES FUNCTIONS




# Create your views here.

#renders the INTRODUCTION page
def intro(request):
    return render(request, 'primsa/intro.html')

#renders the INPUT page
def index(request):
    return render(request, 'primsa/index.html')

#renders the summary of the results from the user's query on Entrez database
def ncbi(request):
    results_esearch={}
    if request.method == 'POST':
        try:
            user_email = request.POST.get('userEmail')
            gene_name= request.POST.get('geneName')
            org_genus= request.POST.get('orgGenus')
            org_kingdom= request.POST.get('orgKingdom')

            #creates a search term based on user's inputs
            search_term= gene_name+'[GENE]' + ' AND ' + org_kingdom+'[filter]' +' AND ' + 'biomol_mrna[PROP]' + ' OR ' + org_genus +'[ALL FIELDS]' + " AND refseq[filter]"

            #Entrez preparation
            Entrez.email=user_email

            #access Nucleotide database using Entrez's esearch method, allows for use of History feature
            handle_esearch=Entrez.esearch('nuccore',search_term, useHistory=True)
            esearch_records=Entrez.read(handle_esearch)

            #parses the esearch_records for Output
            for key, value in esearch_records.items():
                if key == 'Count':
                    results_esearch[key]= value
                elif key=='QueryTranslation':
                    results_esearch[key]=value
                elif key=='WebEnv':
                    results_esearch[key]=value
                elif key=='QueryKey':
                    results_esearch[key]=value

            return render(request, 'primsa/summary_esearch.html', {'data': results_esearch})
        except:
            return HttpResponse ('Something went wrong. Please go back and try again.')

#downloads the sequences from NCBI using Entrez Utils webenv and history features
#allows the user to select the sequences for MSA, and Primer Design
def download(request):
    seq_inf={}
    if request.method=='POST':
        count= request.POST.get('Count')
        count = int(count)
        web_env= request.POST.get('WebEnv')
        query_key= request.POST.get('QueryKey')
        batch_size= 50

        handle_write=open("input_0.fasta",'w')

        #code copied from Biopython tutorial pdf
        for start in range(0, count, batch_size):
                   end=min(count, start+batch_size)
                   #print("Going to download record %i to %i" % (start+1, end))
                   fetch_handle= Entrez.efetch(db='nuccore', rettype="fasta", retstart=start, retmax= batch_size, webenv=web_env, query_key=query_key)
                   data=fetch_handle.read()# can't use Entrez.read(fetch_handle) because the file is not in XML
                   fetch_handle.close()
                   handle_write.write(data)
        handle_write.close()

        handle_open=open("input_0.fasta", 'r')

        for seq_record in SeqIO.parse(handle_open,'fasta'):
            temp=[]
            temp.append(seq_record.seq[0:20])
            temp.append(len(seq_record.seq))
            temp.append( seq_record.description)
            seq_inf[seq_record.id]= temp

        return render(request, 'primsa/retrieve.html', {'data': seq_inf})





#this works: tested April 6th, 2017
def clustal(request):
    #sequence IDs selected by user for MSA
    selected_ids=[]

    #Further screening to remove sequences with size 0 or any genome sequences
    filtered_ids=[]



    #the POST request contains the ids for the selected sequences
    if request.method=='POST':
        for key in request.POST:
            if key!='csrfmiddlewaretoken':
                selected_ids.append(key)

    #FILE HANDLES for reading and writing files
    handle_open= open('input_0.fasta','r')
    handle_write=open('input.fasta','w')

    for seq_record in SeqIO.parse(handle_open,'fasta'):
        #CONSTRAINTS: Sequence Length >0, Sequence NOT Genomic, Sequence SELECTED by the USER
        if (len(seq_record.seq) != 0) and not ('genome' in seq_record.description) and (seq_record.id in selected_ids):
            #CREATES a NEW Fasta file with only the sequences matching the
            #CONSTRAINTS mentioned above
            SeqIO.write(seq_record, handle_write,"fasta")
            filtered_ids.append(seq_record.id)
            split_string= seq_record.description.split()
            first_four= " ".join(split_string[0:4])
            tree_map[seq_record.id]= first_four

    handle_write.close()
    handle_open.close()
    #summary['selected']= selected_ids
    #summary['filtered']= filtered_ids

    results=tree_draw_newick()
    return render(request,'primsa/results.html', {'data':results})

    #testing for changing the node label
    #return render(request,'primsa/summary_esearch.html', {'data':tree_map})










def msa(request):
    #seqID={}
    #seqWithID={}

    '''
    randBaseFileName= ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(15))
    fastaFileName= randBaseFileName+ '.fasta'
    alignFileName= randBaseFileName + '.aln'
    treeFileName=randBaseFileName + '.dnd'
    treeImageName= randBaseFileName + '.jpg'
    basePathImages= '../primsa/static/images/'
    treeImagePath= basePathImages + treeImageName
    '''

    if request.method=='POST':
        write_handle=open("input.fasta",'w')
        write_handle.write(request.POST.get('fastaSeqs'))
        write_handle.close()

        handle_open=open("input.fasta",'r')
        for seq_record in SeqIO.parse(handle_open,"fasta"):
            split_string= seq_record.description.split()
            first_four= " ".join(split_string[0:4])
            tree_map[seq_record.id]= first_four

        results=tree_draw_newick()

        return render(request,'primsa/results.html', {'data':results})



'''
multipleSeqAlignment('filtered.fasta')
alignment= AlignIO.read('filtered.aln', "clustal")
numberCol=alignment.get_alignment_length()
numberRow=len(alignment)
consensusSeq=''
for each in list(range(numberCol)):
    column=alignment[:,each]

    consensusSeq=consensusSeq + countIndBase(column, numberRow)

seq_args={'SEQUENCE_ID':'CONSENSUS_SEQ', 'SEQUENCE_TEMPLATE': consensusSeq}
global_args= {'PRIMER_OPT_SIZE': 20,'PRIMER_PICK_INTERNAL_OLIGO': 1,'PRIMER_PICK_LEFT_PRIMER':1,'PRIMER_PICK_RIGHT_PRIMER':1,'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],}


results=bindings.designPrimers(seq_args, global_args)

return render(request,'primsa/selected_vs_filtered.html', {'data':results})
'''
