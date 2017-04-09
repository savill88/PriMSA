from django.shortcuts import render, HttpResponse
import string, random
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
import pylab
from primer3 import bindings


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
    cline=ClustalwCommandline("clustalw", infile=fileName, output="fasta")
    cline()

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


def download(request):
    seq_inf={}
    if request.method=='POST':
        count= request.POST.get('Count')
        count = int(count)
        web_env= request.POST.get('WebEnv')
        query_key= request.POST.get('QueryKey')
        batch_size= 50
        file_name= 'test.fasta'
        handle_write=open(file_name,'w')

        #code copied from Biopython tutorial pdf
        for start in range(0, count, batch_size):
                   end=min(count, start+batch_size)
                   #print("Going to download record %i to %i" % (start+1, end))
                   fetch_handle= Entrez.efetch(db='nuccore', rettype="fasta", retstart=start, retmax= batch_size, webenv=web_env, query_key=query_key)
                   data=fetch_handle.read()# can't use Entrez.read(fetch_handle) because the file is not in XML
                   fetch_handle.close()
                   handle_write.write(data)
        handle_write.close()

        handle_open=open(file_name, 'r')

        for seq_record in SeqIO.parse(handle_open,'fasta'):
            temp=[]
            temp.append(seq_record.seq[0:20])
            temp.append(len(seq_record.seq))
            temp.append( seq_record.description)
            seq_inf[seq_record.id]= temp

        return render(request, 'primsa/retrieve.html', {'data': seq_inf})





#this works: tested April 6th, 2017
def clustal(request):
    selected_ids=[]
    filtered_ids=[]
    summary={}

    #the POST request contains the ids for the selected sequences
    if request.method=='POST':
        for key in request.POST:
            if key!='csrfmiddlewaretoken':
                selected_ids.append(key)

    #FILE HANDLES for reading and writing files
    handle_open= open('test.fasta','r')
    handle_write=open('filtered.fasta','w')

    for seq_record in SeqIO.parse(handle_open,'fasta'):
        #CONSTRAINTS: Sequence Length >0, Sequence NOT Genomic, Sequence SELECTED by the USER
        if (len(seq_record.seq) != 0) and not ('genome' in seq_record.description) and (seq_record.id in selected_ids):
            #CREATES a NEW Fasta file with only the sequences matching the
            #CONSTRAINTS mentioned above
            SeqIO.write(seq_record, handle_write,"fasta")
            filtered_ids.append(seq_record.id)

    handle_write.close()
    handle_open.close()
    summary['selected']= selected_ids
    summary['filtered']= filtered_ids

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













def msa(request):
    seqID={}
    seqWithID={}
    randBaseFileName= ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(15))
    fastaFileName= randBaseFileName+ '.fasta'
    alignFileName= randBaseFileName + '.aln'
    treeFileName=randBaseFileName + '.dnd'
    treeImageName= randBaseFileName + '.jpg'
    basePathImages= '../primsa/static/images/'
    treeImagePath= basePathImages + treeImageName
    #return HttpResponse(treeFileName)
    if request.method=='POST':
        write_handle=open(fastaFileName,'w')
        write_handle.write(request.POST.get('fastaSeqs'))
        write_handle.close()
        multipleSeqAlignment(fastaFileName)
        #open_handle=open(fastaFileName,'r')
        alignment= AlignIO.read(alignFileName, "clustal")
        numberCol=alignment.get_alignment_length()
        numberRow=len(alignment)


        consensusSeq=''
        for each in list(range(numberCol)):
            column=alignment[:,each]

            consensusSeq=consensusSeq + countIndBase(column, numberRow)
        #return HttpResponse(consensusSeq)
        seq_args={'SEQUENCE_ID':'CONSENSUS_SEQ', 'SEQUENCE_TEMPLATE': consensusSeq}
        global_args= {'PRIMER_OPT_SIZE': 20,'PRIMER_PICK_INTERNAL_OLIGO': 1,'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],}


        results=bindings.designPrimers(seq_args, global_args)

        treeNewick= Phylo.read(treeFileName, "newick")
        #Phylo.draw_ascii(treeNewick)
        treeNewick.rooted=True
        #note sure what this does
        treeNewick.ladderize()
        Phylo.draw(treeNewick, do_show=False)
        pylab.axis('off')
        pylab.savefig("/home/savill88/Desktop/senior-project/PriMSA/django_PriMSA/primsa/static/images/tree", format='png', bbox_inches='tight', dpi=300)


        #results= bindings.runP3Design()

        return render(request,'primsa/selected_vs_filtered.html', {'data':results})

        '''
        treeNewick= Phylo.read(treeFileName, "newick")
        #Phylo.draw_ascii(treeNewick)
        treeNewick.rooted=True
        #note sure what this does
        treeNewick.ladderize()
        Phylo.draw(treeNewick, do_show=False)
        pylab.axis('off')
        pylab.savefig("/home/savill88/Desktop/senior-project/PriMSA/django_PriMSA/primsa/static/images/tree", format='png', bbox_inches='tight', dpi=300)

        return HttpResponse('Whatever')
        '''
        '''
        for seq_record in SeqIO.parse(open_handle,"fasta"):
            temp=[]
            temp.append(seq_record.seq[0:20])
            temp.append(len(seq_record.seq))
            seqID[seq_record.id]= temp
        return render(request, 'primsa/retrieve.html', {'data': seqID})
        #return HttpResponse(fastaInput)
        '''
        #return HttpResponse('Did not enter the POST if ')
