from django.shortcuts import render, HttpResponse
import requests
import json
from Bio import Entrez

# Create your views here.

def index(request):
    #return HttpResponse('This is the index page!!')
    return render(request, 'primsa/index.html')
