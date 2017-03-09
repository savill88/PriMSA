from django.shortcuts import render, HttpResponse
import requests
import json
from Bio import Entrez, SeqIO

# Create your views here.

def index(request):

    return render(request, 'primsa/index.html', {})

    #return render(request, 'primsa/index.html')



    #return HttpResponse('This is the index page!!')
    #return render(request, 'primsa/index.html')


def retrieve(request):
    parsedData=[]

    if request.method == 'POST':
        #saves the user's input
        userEmail = request.POST.get('userEmail')
        geneName= request.POST.get('geneName')
        orgGenus= request.POST.get('orgGenus')
        orgKingdom= request.POST.get('orgKingdom')

        #creates a search term based on user's inputs
        searchTerm= geneName+'[GENE]' + ' AND ' + orgKingdom+'[ORG]' +' AND ' + 'biomol_mrna[PROP]'

        #Entrez preparation
        Entrez.email=userEmail

        #access Nucleotide database using Entrez's esearch method
        handle=Entrez.esearch('nuccore',searchTerm)
        record=Entrez.read(handle)
        parsedData.append(record['IdList'])


    return render(request, 'primsa/search.html', {'data': parsedData})

def msa(request):
    seqID={}
    if request.method=='POST':
        write_handle=open('test.fasta','w')
        write_handle.write(request.POST.get('fastaSeqs'))
        write_handle.close()
        open_handle=open('test.fasta','r')
        for seq_record in SeqIO.parse(open_handle,"fasta"):
            temp=[]
            temp.append(seq_record.seq[0:20])
            temp.append(len(seq_record.seq))
            seqID[seq_record.id]= temp
        return render(request, 'primsa/retrieve.html', {'data': seqID})
        #return HttpResponse(fastaInput)

    return HttpResponse('Did not enter the POST if ')
