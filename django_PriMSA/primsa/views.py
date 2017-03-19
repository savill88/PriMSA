from django.shortcuts import render, HttpResponse
import string, random
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
import pylab
from primer3 import bindings


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
    try:
        if request.method == 'POST':
            user_email = request.POST.get('userEmail')
            gene_name= request.POST.get('geneName')
            org_genus= request.POST.get('orgGenus')
            org_kingdom= request.POST.get('orgKingdom')

            #creates a search term based on user's inputs
            search_term= gene_name+'[GENE]' + ' AND ' + org_kingdom+'[ORG]' +' AND ' + 'biomol_mrna[PROP]' + ' OR ' + org_genus +'[ALL FIELDS]'

            #Entrez preparation
            Entrez.email=user_email

            #access Nucleotide database using Entrez's esearch method, allows for use of History feature
            handle_esearch=Entrez.esearch('nuccore',search_term, useHistory=True)
            esearch_records=Entrez.read(handle_esearch)

            #parses the esearch_records for Output
            for key, value in esearch_records.items():
                if key == 'Count':
                    results_esearch[key]= value
                elif key=='IdList':
                    results_esearch[key]=value
                elif key=='QueryTranslation':
                    results_esearch[key]=value

            return render(request, 'primsa/summary_esearch.html', {'data': results_esearch})
    except:
        print "Could not perform the search!!"

def countIndBase(col, numRows):
    count = {'A' : 0, 'T' : 0, 'C' : 0, 'G' : 0, 'N' : 0, '-':0, 'W':0,'S':0,'H':0, 'Y':0,'M':0,'R':0}
    for letter in col:
        count[letter] += 1

    #print count
    '''
    for key in count:
    	print count[key]
    	#= count[key]/numRows

    print count
    '''
    baseNuc= max(count, key=count.get)
    if (baseNuc=='-'):
        return 'N'
    #print " The returning base is: ", baseNuc
    else:
        return max(count, key=count.get)









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
        #open_handle=open(fastaFileName,'r')
        cline=ClustalwCommandline("clustalw", infile=fastaFileName)
        cline()
        alignment= AlignIO.read(alignFileName, "clustal")
        numberCol=alignment.get_alignment_length()
        numberRow=len(alignment)


        consensusSeq=''
        for each in list(range(numberCol)):
            column=alignment[:,each]
            #print each, "Col Number"
            consensusSeq=consensusSeq + countIndBase(column, numberRow)
        #return HttpResponse(consensusSeq)
        seq_args={'SEQUENCE_ID':'CONSENSUS_SEQ', 'SEQUENCE_TEMPLATE': consensusSeq}
        global_args= {'PRIMER_OPT_SIZE': 20,'PRIMER_PICK_INTERNAL_OLIGO': 1,'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],}


        results=bindings.designPrimers(seq_args, global_args)

        #results= bindings.runP3Design()

        return HttpResponse(results)

        '''
        treeNewick= Phylo.read(treeFileName, "newick")
        #Phylo.draw_ascii(treeNewick)
        treeNewick.rooted=True
        #note sure what this does
        treeNewick.ladderize()
        Phylo.draw(treeNewick, do_show=False)
        pylab.axis('off')
        pylab.savefig(treeImagePath, format='jpg', bbox_inches='tight', dpi=300)

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
